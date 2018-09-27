#include "IHSOptimizer.h"
#include "CPLEXInterface.h"

#include <unordered_set>
#include <vector>

#include "WaspFacade.h"
#include "outputBuilders/MultiOutputBuilder.h"

using namespace std;

IHSOptimizer::IHSOptimizer(ostream & log)
    : log(log),
      ASP_solver(nullptr),
      IP_solver(nullptr),
      avgCoreSize(0), avgCoreSizeMinimized(0),
      LB(0), UB(numeric_limits<uint64_t>::max()),
      cfg(GlobalConfig::get())
{
}

IHSOptimizer::~IHSOptimizer() {
  delete IP_solver;
  delete ASP_solver;
}

void IHSOptimizer::solve(istream & instance) {

  solveTimer.start();

  ASP_solver = new WaspFacade();
  IP_solver = new CPLEXInterface(log);

  ASP_solver->attachProgramListener(this);

  ASP_solver->setOutputPolicy(MULTI);
  /*
  if (cfg.silent)
    ASP_solver->disableOutput();
  */
  ASP_solver->setRestartsPolicy(SEQUENCE_BASED_RESTARTS_POLICY, 100);

  ASP_solver->readInput(instance);
  //ASP_solver->setMinimizeUnsatCore(cfg.minimizeCores);

  unsigned n_vars = ASP_solver->numberOfVariables();
  unsigned n_levels = ASP_solver->numberOfLevels();

  log << "c #Vars " << n_vars << endl;
  log << "c #Levels " << n_levels << endl;

  vector< uint64_t > levelWeights;

  //
  // Gather set of all optimization literals
  // and compute total weight for each level
  //
  for (unsigned level = 0; level < n_levels; ++level) {
    vector< OptimizationLiteralData > tmp;
    ASP_solver->getOptLiterals(level, tmp);
    optLits.insert(optLits.end(), tmp.begin(), tmp.end());

    uint64_t levelWeight = 0;

    for (auto data : tmp) {
      uint64_t litWeight = data.weight;
      if (level > 0) {
        litWeight += levelWeights[level - 1];
      }
      levelWeight += litWeight;
    }

    levelWeights.push_back(levelWeight);
  }

  //
  // Add optimization literals to IP solver
  // Modify weight if necessary due to level
  //
  for (auto & data : optLits) {
    uint64_t weight = data.weight;
    if (data.level > 0) {
      weight += levelWeights[data.level - 1];
    }

    log << "c Literal " << data.lit << "\t"
        << "Weight " << data.weight << "\t"
        << "AdjWeight " << weight << "\t"
        << "Level " << data.level << "\t"
        << "ID " << data.lit.getId() << endl;

    data.weight = weight;

    id_weight[data.lit.getId()] = weight;
    IP_solver->addObjectiveVariable(data.lit.getId(), weight);
  }


  vector< Literal > core;
  vector< int > coreIds;
  vector< Literal > assumptions;
  unordered_set<int> hittingSet;

  uint64_t init_UB = 0;
  for (auto data : optLits)
    init_UB += data.weight;
  updateUB(init_UB);

  enum hsPhase { optimal, greedy, disjoint };

  hsPhase phase = optimal;
  bool modelExists = false;

  unsigned nonoptCoreCount = 0;

  unordered_set<int> fixedVars;

  if (cfg.clarkeIpLB || cfg.clarkeLpLB) {
    // create temporary IP instance of clarke's completion
    CPLEXInterface *clarke_IP = new CPLEXInterface(log);
    for (auto pair : id_weight)
      clarke_IP->addObjectiveVariable(pair.first, pair.second);

    feedClauses(clarke_IP);
 
    // solve IP or LP relaxation
    unordered_map<int, bool> forcedVars;
    if (cfg.clarkeIpLB) {
      clarke_IP->computeHS(hittingSet, LB, UB, forcedVars);
    } else { // clarkeLpLB
      LB = clarke_IP->computeRelaxedBound(LB, UB, forcedVars);
    }

    log << "c clarke LB " << LB << endl;

    delete clarke_IP;
  }

  if (cfg.ipClauses) feedClauses(IP_solver);

  if (cfg.persistLpLB) feedClauses(IP_solver, true);

  if (cfg.ipSeeding) eqSeeding();

  //
  // Loop until an optimal solution is found
  //
  int res;
  for (unsigned iteration = 0;;iteration++) {

    log << endl << "c iteration " << iteration << " ["<<LB<<".."<<UB<<"] ("
        << solveTimer.cpu_ms_total() << " ms)" << endl;

    if (LB > UB) exit(1);

    unordered_map<int, bool> forcedVars;

    switch (phase) {
      case optimal:

        IP_solver->computeHS(hittingSet, LB, UB, forcedVars);

        if (LB == UB && modelExists)
          goto optimum_found;

        for (auto pair : forcedVars) {
          //log << "c HS LB forces " << pair.first << "=" << pair.second << endl;
          vector<Literal> unit = {
            Literal::createLiteralFromInt(pair.first * (pair.second ? 1 : -1))
          };
          fixedVars.insert(abs(pair.first));
          ASP_solver->addClause(unit);
        }

        log << "c Optimal HS";
        log << " weight=" << LB;
        break;
      case greedy:
        log << "c Greedy HS";
        hittingSet = greedyHittingSet();
        break;
      case disjoint:
        log << "c Disjoint HS";
        hittingSet.insert(coreIds.begin(), coreIds.end());
        break;
    }
    log << " size=" << hittingSet.size() << endl;

    setAssumptions(hittingSet, assumptions);

    core.clear();
    coreIds.clear();
    coreExtractionTimer.start();
    res = ASP_solver->solve(assumptions, core);
    coreExtractionTimer.stop();

    if (res == COHERENT) {
      //
      // If hitting set was minimum, we have an optimal solution
      // Otherwise, break disjoint phase loop and
      // compute the minimum hitting set next iteration
      //
      // log << "c ANSWER SET" << endl;

      uint64_t new_cost = answerSetCost();
      if (new_cost < UB || !modelExists) {
        ASP_solver->printAnswerSet();
        ASP_solver->printOptimizationValue();
        modelExists = true;
        updateUB(new_cost);
      }

      if (LB == UB && modelExists)
        goto optimum_found;

      switch (phase) {
        case optimal:
          // optimal solution
          log << "c IP opt " << LB << endl;
          ASP_solver->printAnswerSet();
          goto optimum_found;
        case greedy:
          phase = optimal;
          break;
        case disjoint:
          if (cfg.greedyHittingSets)
            phase = greedy;
          else
            phase = optimal;
          break;
      }

    } else if (res == INCOHERENT) {

      if (core.size() == 0) {
        // core is empty, incoherence not due to assumptions
        ASP_solver->printIncoherence();
        cout << "INCOHERENT" << endl;
        return;
      }

      if (cores.size()) {
        avgCoreSize = (avgCoreSize * cores.size() + core.size()) / (cores.size() + 1);
      } else {
        avgCoreSize = core.size();
      }

      if (cfg.minimizeCores) {
        minimizeCore(core);
      } else if (cfg.cardinalityMinimize) {
        cardMinimize(core);
      }

      // log << "c INCOHERENCE" << endl;
      log << "c core size " << core.size() << endl;

      for (auto lit : core)
        coreIds.push_back(-lit.getId());

      IP_solver->addCore(coreIds);
      cores.push_back(coreIds);

      for (auto id : coreIds) {
        if ( !id_coreCount.count(id) ) {
          id_coreCount[id] = 1;
        } else {
          id_coreCount[id]++;
        }
      }

      switch (phase) {
        case optimal:
          if (cfg.disjointPhases)
            phase = disjoint;
          else if (cfg.greedyHittingSets)
            phase = greedy;
          break;
        case greedy:
          nonoptCoreCount++;
          if (cfg.disjointPhases)
            phase = disjoint;
          break;
        case disjoint:
          nonoptCoreCount++;
          break;
      }

      if (nonoptCoreCount >= cfg.nonoptCoreLimit) {
        nonoptCoreCount = 0;
        phase = optimal;
      }
    }

    continue;
    optimum_found: break;
  } // end main loop

  //ASP_solver->enableOutput();

  ASP_solver->printOptimumFound();

  ASP_solver->getOutputBuilder()->onFinish();
}

void IHSOptimizer::updateUB(uint64_t cost) {
  if (cost < UB) {
    log << "c UB " << cost << endl;
    UB = cost;
  }
}

void IHSOptimizer::setAssumptions(unordered_set< int > & hittingSet, vector< Literal > & assumptions) {
  //
  // If an optimization literal is not in the hitting set
  // add its opposite literal to the list of assumptions,
  // forcing the corresponding constraint to be satisfied
  //
  assumptions.clear();

  //log << "c Assumptions";
  for (auto data : optLits) {
    if ( ! hittingSet.count(data.lit.getId()) ) {
      assumptions.push_back(data.lit.getOppositeLiteral());
      //log << " " << -data.lit.getId();
    }
  }
  //log << endl;
}

uint64_t IHSOptimizer::answerSetCost() {
  //
  // Compute sum of true optimization literals in latest answer set
  //
  uint64_t cost = 0;
  for (auto data : optLits) {
    if (ASP_solver->isTrue(data.lit)) {
      cost += data.weight;
    }
  }

  return cost;
}

void IHSOptimizer::foundCore(vector< Literal > & core) {

  vector< int > coreIds;
  log << "c INCOHERENCE" << endl;
  log << "c core";
  for (auto lit : core) {
    log << " " << lit.getId();
    coreIds.push_back(-lit.getId());
  }
  log << endl;

  IP_solver->addCore(coreIds);
  cores.push_back(coreIds);

  for (auto id : coreIds) {
    if ( !id_coreCount.count(id) ) {
      id_coreCount[id] = 1;
    } else {
      id_coreCount[id]++;
    }
  }
}

void IHSOptimizer::minimizeCore(vector< Literal > & core) {

  vector< Literal > criticalLits;
  vector< Literal > untestedLits(core);

  vector< Literal > assumptions;

  minimizeTimer.start();
  bool minimal = true;

  while (untestedLits.size()) {

    Literal testLit = untestedLits.back();
    untestedLits.pop_back();

    assumptions.clear();

    //log << "c Assumptions";
    //for (auto lit : criticalLits)
    //  log << " " << lit.getId();
    //for (auto lit : untestedLits)
    //  assumptions.push_back(lit);
    //log << endl;

    for (auto lit : criticalLits)
      assumptions.push_back(lit);

    for (auto lit : untestedLits)
      assumptions.push_back(lit);

    vector< Literal > subCore;

    ASP_solver->setBudget(BUDGET_TIME, cfg.minimizeCallSeconds);

    int res = ASP_solver->solve(assumptions, subCore);

    if (res == COHERENT) {
      //log << "Core sat: Lit critical" << endl;
      criticalLits.push_back(testLit);

    } else if (res == INCOHERENT) {

      //log << "Core unsat: Lit redundant" << endl;
      //log << "Subcore size " << subCore.size() << endl;
      //log << "Assumptions size " << assumptions.size() << endl;

      //  discard testLit
    } else if (res == INTERRUPTED) {
      cout << "int" << endl;
      criticalLits.push_back(testLit);
      minimal = false;
    }

    if (minimizeTimer.cpu_ms_total() > cfg.minimizeTimeLimit * solveTimer.cpu_ms_total()) {
      for (auto lit : untestedLits) {
        criticalLits.push_back(lit);
      }
      minimal = false;
      break;
    }
  }

  log << "c Minimized: " << core.size() << " -> " << criticalLits.size();
  if (minimal) cout << endl;
  else cout << " (*)" << endl;

  core.swap(criticalLits);

  if (cores.size()) {
    avgCoreSizeMinimized = (avgCoreSizeMinimized * cores.size() + core.size()) / (cores.size() + 1);
  } else {
    avgCoreSizeMinimized = core.size();
  }

  ASP_solver->budgetOff();

  minimizeTimer.stop();
}

Literal IHSOptimizer::addTemporaryCard(vector< Literal > & core, unsigned int bound) {

  vector< Literal > lits(core);

  // add fresh variable
  Var id_var = ASP_solver->addVariable();
  Literal aux(id_var);
  lits.push_back( aux );

  vector< uint64_t > weights;
  for (unsigned i = 0; i < core.size(); ++i) weights.push_back(1);

  // give new variable bound weight
  weights.push_back(bound);

  // constraint trivially satisfied by aux = 1
  ASP_solver->addPseudoBooleanConstraint(lits, weights, bound);

  return aux;
}

void IHSOptimizer::cardMinimize(vector< Literal > & core) {

  unsigned int oldSize = core.size();

  minimizeTimer.start();

  vector< Literal > mus;

  //cout << "core size " << oldSize << endl;

  while (true) {

    Literal aux = addTemporaryCard(core, core.size() - 1);
    vector< Literal > assumptions({ aux.getOppositeLiteral() });

    vector< Literal > subCore;

    for (auto lit : mus)
      assumptions.push_back(lit);

    int res = ASP_solver->solve(assumptions, subCore);

    if (res == COHERENT) {
      //log << core.size() - 1 << " satisfiable" << endl;

      for (unsigned int i = 0; i < core.size(); ++i) {
        if (ASP_solver->isFalse(core[i])) {
          //log << core[i].getId() << " was false"  << endl;
          mus.push_back(core[i]);
          core[i] = core[core.size() - 1];
          core.pop_back();
          break;
        }
      }
    } else if (res == INCOHERENT) {

      if (find(subCore.begin(), subCore.end(), aux.getOppositeLiteral()) != subCore.end()
          && core.size()) {
        //log << "card in core" << endl;
        core.pop_back();
      } else {
        //log << "mus is unsat" << endl;
        ASP_solver->addClause( aux );
        core.swap(mus);
        break;
      }
      //log << "Core sat: Lit critical" << endl;
      //log << core.size() - 1 << " unsatisfiable" << endl;
      //core.swap(subCore);
    }

    // deactivate card constraint
    ASP_solver->addClause( aux );
  }

  log << "c Minimized: " << oldSize << " -> " << core.size() << endl;

  if (cores.size()) {
    avgCoreSizeMinimized = (avgCoreSizeMinimized * cores.size() + core.size()) / (cores.size() + 1);
  } else {
    avgCoreSizeMinimized = core.size();
  }

  minimizeTimer.stop();

}

/*
void IHSOptimizer::minimizeCore(vector< Literal > & core) {

  vector< Literal > assumptions;

  minimizeTimer.start();

  int a = 0;
  int b = core.size() - 1;

  vector< Literal > subCore;

  //cout << "CORE " << core.size() << endl;

  // find smallest j for which core[0..j] is unsat
  while (a < b) {
    //cout << a << " " << b << endl;
    int k = (a + b) / 2;

    assumptions.clear();

    // test literals 0 .. j+c
    for (unsigned i = 0; i < k; ++i) {
      assumptions.push_back(core[i]);
    }

    int res = ASP_solver->solve(assumptions, subCore);

    if (res == COHERENT) {
      //cout << "COHERENT " << assumptions.size() << endl;
      a = k+1;
    } else { //if (res == INCOHERENT) {
      //cout << "INCOHERENT " << assumptions.size() << endl;
      b = k;
    }
  }

  //cout << "END " << assumptions.size() << endl;

  log << "c Minimized: " << core.size() << " -> " << b+1 << endl;

  while (core.size()-1 > b) {
    core.pop_back();
  }

  //cout << "MINIMIZED " << core.size() << endl;

  if (cores.size()) {
    avgCoreSizeMinimized = (avgCoreSizeMinimized * cores.size() + core.size()) / (cores.size() + 1);
  } else {
    avgCoreSizeMinimized = core.size();
  }

  minimizeTimer.stop();

}
*/

void IHSOptimizer::printStats() {
  log << "c cores found: " << cores.size() << endl;
  log << "c average core size: " << avgCoreSize << endl;
  if (cfg.minimizeCores || cfg.cardinalityMinimize)
    log << "c minimized core size: " << avgCoreSizeMinimized << endl;
  log << "c core extraction time: " << coreExtractionTimer.cpu_ms_total() << " ms" << endl;
  if (cfg.minimizeCores)
    log << "c minimization time: " << minimizeTimer.cpu_ms_total() << " ms" << endl;

  IP_solver->printStats();
}

unordered_set< int > IHSOptimizer::greedyHittingSet() {
  unordered_set< int > hs;

  // keep track of the number of occurrences of each variable in the cores
  // and also the variable weight, to reduce lookups
  unordered_map<int, pair<int, uint64_t>> var_count_weight;

  // remember which cores each variable occurs in
  unordered_map<int, vector<int>> var_cores;

  log << "c finding greedy hitting set of " << cores.size() << " cores" << endl;

  for (unsigned i = 0; i < cores.size(); ++i) {
    for (int v : cores[i]) {
      if (var_cores.find(v) == var_cores.end()) {
        var_cores[v] = vector<int>();
        var_count_weight[v] = make_pair(id_coreCount[v], id_weight[v]);
      }
      var_cores[v].push_back(i);
    }
  }

  // keep track of which cores are hit by the current vector hs
  vector<bool> hit(cores.size(), false);

  unsigned unhits = cores.size();
  uint64_t weight = 0;

  // repeat until all cores are hit
  while (unhits) {

    // find the variable v with lowest (weight / #occurrences)
    int v = -1;
    double o = std::numeric_limits<double>::max();
    uint64_t w = std::numeric_limits<uint64_t>::max();

    for (auto& e : var_count_weight) {
      int count = e.second.first;
      uint64_t wt = e.second.second;
      double _o = ((double)wt / (double)count);
      if (count != 0 && _o < o) {
        o = _o;
        w = wt;
        v = e.first;
      }
    }

    weight += w;

    // update variable counts
    for (int c : var_cores[v]) {
      // if that core hasn't been hit yet
      if (!hit[c]) {
        // for each variable in the core, decrement occurrences
        for (int b : cores[c]) {
          var_count_weight[b].first -= 1;
        }
      }
    }

    // track hit cores
    for (int c : var_cores[v]) {
      if (!hit[c]) {
        hit[c] = true;
        --unhits;
      }
    }

    hs.insert(v);
  }

  return hs;
}

//
// ProgramListener methods
//

void IHSOptimizer::addedClause(const vector<int>& literals) {
  /*
  log << "c listener clause";
  for (int l : literals) {
    log << " " << l;
  }
  log << endl;
  */
  clauses.push_back(literals);
}

void IHSOptimizer::addedAggregate(unsigned int id,
  const vector<int>& literals,
  const vector<uint64_t>& weights,
  uint64_t bound) {
  /*
  log << "c listener aggregate " << id << ": ";
  for (unsigned i = 0; i < literals.size(); ++i) {
    log << " " << literals[i] << "(" << weights[i] << ")";
  }
  log << " <= " << bound << endl;
  */
  Aggregate agg = {id, literals, weights, bound};

  aggregates.push_back(agg);
}

void IHSOptimizer::eqSeeding() {
  //
  // Eq Seeding etc
  //
  log << "c #Aggregates: " << aggregates.size() << endl;
  log << "c #Clauses: " << clauses.size() << endl;

  unsigned aggregatesWithSoftLiteral = 0;

  unordered_set<int> softVariables;
  unordered_map<int, int> softClauseCount;
  unordered_map<int, int> opt_eq_plain;

  // Create set of optimization literals and initialize clause counters
  for (auto data : optLits) {
    softVariables.insert(data.lit.getVariable());
    softClauseCount[data.lit.getVariable()] = 0;
  }

  // Find optimization literals appearing in exactly one clause
  for (auto cl : clauses) {

    for (unsigned i = 0; i < cl.size(); ++i) {
      int l = cl[i];  // opt lit
      bool s = l < 0; // sign
      int v = abs(l); // opt var

      if (softVariables.find(v) != softVariables.end()) {

        // huom: post-increment
        unsigned count = softClauseCount[v]++;

        // variable was seen in another clause, is not an equivalence
        if (count) opt_eq_plain.erase(v);

        // binary clause with new optimization literal?
        if (cl.size() == 2 && !count) {
          int l_other = cl[cl.size() - 1 - i];
          opt_eq_plain[v] = s ? l_other : -l_other;

          /*
          cout << "c eq from cl";
          for (auto c : cl) cout << " " << c;
          cout << endl;
          */
        }
      }
    }
  }

  // Check aggregates, remove equivalence if soft literal present
  for (auto agg : aggregates) {
    for (int l : agg.literals) {
      int v = abs(l);
      if (softVariables.find(v) != softVariables.end()) {
        aggregatesWithSoftLiteral++;
        opt_eq_plain.erase(v);
      }
    }
  }

  unordered_map<int,int> plain_eq_opt;

  // reverse opt <-> plain map and add negations
  for (auto opt_plain : opt_eq_plain) {
    plain_eq_opt[opt_plain.second] = opt_plain.first;
    plain_eq_opt[-opt_plain.second] = -opt_plain.first;
  }

  /*
  for (auto p : plain_eq_opt) {
    cout << p.first << " <-> " << p.second << endl;
  }
  */

  // add also self-equivalences for opt lits
  for (int v : softVariables) {
    plain_eq_opt[v] = v;
    plain_eq_opt[-v] = -v;
  }


  // check all clauses:
  // if each literal found in plain<->opt
  // then clause can be seeded as IP constraint.
  vector<vector<int>> eqConstraints;

  for (auto cl : clauses) {
    vector<int> eqConstr;
    for (int l : cl) {
      if (plain_eq_opt.find(l) != plain_eq_opt.end()) {
        eqConstr.push_back(plain_eq_opt[l]);
       } else {
        goto skip_cl;
      }
    }
    eqConstraints.push_back(eqConstr);
    /*
    cout << "c cl";
    for (auto c : cl) cout << " " << c;
    cout << endl;

    cout << "c eq";
    for (auto c : eqConstr) cout << " " << c;
    cout << endl;
    */
    skip_cl: continue;
  }

  unsigned eq_aggs = 0;

  for (auto agg : aggregates) {

    if (plain_eq_opt.find(agg.id) == plain_eq_opt.end())
      goto skip_agg;

    for (int l : agg.literals) {
      if (plain_eq_opt.find(l) != plain_eq_opt.end()) {
        // ok
       } else {
        goto skip_agg;
      }
    }

    eq_aggs++;

    skip_agg: continue;
  }

  log << "c Eq literals: " << opt_eq_plain.size() << endl;
  log << "c Eq constraints: " << eqConstraints.size() << endl;
  log << "c Eq aggregates: " << eq_aggs << endl;

  for (auto constr : eqConstraints)
    IP_solver->addCore(constr);
}

void IHSOptimizer::feedClauses(CPLEXInterface * solver, bool lp_only) {
  log << "c add " << clauses.size() << " clauses to IP" << endl;

  int count = 0;

  for (auto cl : clauses) {
    solver->addClause(cl, lp_only);
    if ( ++count >= cfg.lpClauseLimit )
      break;
  }

  log << "c add " << aggregates.size() << " aggregates to IP" << endl;
  for (auto agg : aggregates) {
    /*
    log << agg.id << " <->";
    for (int i = 0; i < agg.literals.size(); ++i) {
      log << " " << agg.weights[i] << "*" << agg.literals[i];
      if (i != agg.literals.size() - 1) log << " +";
    }
    log << " >= " << agg.bound << endl;
    */
    solver->addAggregate(agg.id, agg.literals, agg.weights, agg.bound, lp_only);
    count += 2;
    if ( count >= cfg.lpClauseLimit )
      break;
  }
}