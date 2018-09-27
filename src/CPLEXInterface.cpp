#include "CPLEXInterface.h"
#include "IHSOptimizer.h"
#include "Timer.h"

#include <iostream>

#define EPS 1e-7

ILOSTLBEGIN

ILOMIPINFOCALLBACK1(LBCutoffCallback,
                    IloNum, lb)
{
  if ( hasIncumbent() && getIncumbentObjValue() == lb) {
    //log << "c IP LB cutoff " << getIncumbentObjValue() << endl;
    abort();
  }
}

CPLEXInterface::CPLEXInterface(ostream & log) : 
    log(log), 
    objFuncAttached(false), 
    nObjVars(0), 
    nVars(0),
    ipCalls(0), 
    reducedCostForcedVars(0), 
    reducedCostRelaxedVars(0),
    cfg(GlobalConfig::get())  {

  env = new IloEnv();
  model = IloModel(*env);
  objVars = IloNumVarArray(*env);
  vars = IloNumVarArray(*env);
  cons = IloRangeArray(*env);
  cplex = IloCplex(model);

  lp_model = IloModel(*env);
  lp_model.add(model);
  lp_cplex = IloCplex(lp_model);
  lp_model.add(IloConversion(*env, objVars, ILOFLOAT));
  lp_model.add(IloConversion(*env, vars, ILOFLOAT));
  lp_cplex.setParam(IloCplex::Threads, 1);
  lp_cplex.setParam(IloCplex::SimDisplay, 0);
  lp_cplex.setParam(IloCplex::MIPDisplay, 0);

  cplex.setParam(IloCplex::Threads, 1);
  cplex.setParam(IloCplex::EpGap, 0.0);
  cplex.setParam(IloCplex::EpAGap, 0.0);
  cplex.setParam(IloCplex::MIPDisplay, 0);

  objective = IloMinimize(*env, IloExpr(*env));
}

CPLEXInterface::~CPLEXInterface() {
  (*env).end();
  delete env;
}

// Create a CPLEX variable for the bvar and append
// it to the objective function with the given weight
void CPLEXInterface::addObjectiveVariable(int id, uint64_t weight) {
  // detact objective function to modify it
  if (objFuncAttached) {
    model.remove(objective);
    objFuncAttached = false;
  }

  // print variable id to string
  char name[11];
  snprintf(name, 11, "%10d", id);

  // create CPLEX variable
  IloNumVar x(*env, 0, 1, ILOBOOL, name);
  var_to_IloVar[id] = x;
  var_to_weight[id] = weight;
  objVars.add(x);
  nObjVars++;

  //lp_model.add(IloConversion(*env, x, ILOFLOAT));

  // update objective function
  objective.setExpr(objective.getExpr() + double(weight) * x);
}

void CPLEXInterface::addVariable(int id) {
  char name[11];
  snprintf(name, 11, "%10d", id);

  IloNumVar x(*env, 0, 1, ILOBOOL, name);

  var_to_IloVar[id] = x;
  vars.add(x);

  //lp_model.add(IloConversion(*env, x, ILOFLOAT));

  nVars++;
}

void CPLEXInterface::addCore(vector<int>& core) {

  IloExpr expr(*env);

  int neg = 0;
  for (int c : core) {
    if (var_to_IloVar.find(c) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[c];
      expr += x;
    } else if (var_to_IloVar.find(-c) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[-c];
      expr -= x;
      ++neg;
    } else {
      printf("Variable %d not found\n", c);
      exit(1);
    }
  }

  IloRange con = expr >= (1 - neg);

  cons.add(con);
  model.add(con);
}

void CPLEXInterface::addClause(vector<int>& cl, bool lp_only) {
 IloExpr expr(*env);

  int neg = 0;
  for (int l : cl) {
    int v = abs(l);

    if (var_to_IloVar.find(l) == var_to_IloVar.end() and
        var_to_IloVar.find(-l) == var_to_IloVar.end()) {
      addVariable(v);
    }

    if (var_to_IloVar.find(l) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[l];
      expr += x;
    } else if (var_to_IloVar.find(-l) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[-l];
      expr -= x;
      ++neg;
    } else {
      printf("Variable %d not found\n", v);
      exit(1);
    }
  }

  IloRange con = expr >= IloInt(1 - neg);

  cons.add(con);
  if (lp_only)
    lp_model.add(con);
  else
    model.add(con);
}

void CPLEXInterface::addAggregate(int id, std::vector<int> literals, std::vector<uint64_t> weights, uint64_t bound, bool lp_only) {

  IloExpr sum_expr(*env);

  IloInt sum_bound = bound;

  for (unsigned i = 0; i < weights.size(); ++i) {
    int l = literals[i];
    IloInt w = weights[i];

    if (var_to_IloVar.find(l) == var_to_IloVar.end() and
        var_to_IloVar.find(-l) == var_to_IloVar.end()) {
      addVariable(l);
    }

    if (var_to_IloVar.find(l) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[l];
      sum_expr += w*x;
    } else if (var_to_IloVar.find(-l) != var_to_IloVar.end()) {
      IloNumVar x = var_to_IloVar[-l];
      sum_expr -= w*x;

      sum_bound -= w;

    } else {
      printf("Variable %d not found\n", l);
      exit(1);
    }
  }

  if (var_to_IloVar.find(id) == var_to_IloVar.end() and
      var_to_IloVar.find(-id) == var_to_IloVar.end()) {
    addVariable(id);
  }

  IloInt M = sum_bound;
  for (uint64_t w : weights) M += w;

  M *= 2;

  IloNumVar var_id;
  IloRange eq_if;
  IloRange eq_only_if;

  if (var_to_IloVar.find(id) != var_to_IloVar.end()) {
    sum_expr -= M * var_to_IloVar[id];
    eq_if = sum_expr >= IloInt(sum_bound) - M;
    eq_only_if = IloNum(sum_bound - EPS) >= sum_expr;
  } else {
    sum_expr += M * var_to_IloVar[-id];
    eq_if = sum_expr >= IloInt(sum_bound);
    eq_only_if = IloNum(sum_bound + M - EPS) >= sum_expr;
  }

  cons.add(eq_if);
  cons.add(eq_only_if);

  if (lp_only) {
    lp_model.add(eq_if);
    lp_model.add(eq_only_if);
  } else {
    model.add(eq_if);
    model.add(eq_only_if);
  }
}

// returns positive valued objective function variables in 'solution'
// and value of objective function in 'cost'
bool CPLEXInterface::computeHS(unordered_set<int>& out_hittingSet,
    uint64_t & LB, uint64_t UB, unordered_map<int,bool> &forcedVars)
{
  try {

    // attach objective function if detached
    if (!objFuncAttached) {
      model.add(objective);
      objFuncAttached = true;
    }

    out_hittingSet.clear();

    static uint64_t old_LB = 0;
    static uint64_t old_UB = 0;

    bool bounds_changed = old_LB != LB || old_UB != UB;

    // Bounding is enabled and bounds have changed
    if (cfg.reducedCosts && bounds_changed) {
      uint64_t lp_LB = computeRelaxedBound(LB, UB, forcedVars);
      if (lp_LB > LB) LB = lp_LB;
    }

    // stop search if we have previously proved lower bound at current incumbent solution cost
    cplex.use(LBCutoffCallback(*env, LB));

    ++ipCalls;
    ipTimer.start();
    bool ret = cplex.solve();
    ipTimer.stop();

    if (!ret) return false;

    // get objective variable values from cplex
    IloNumArray vals(*env);
    cplex.getValues(vals, objVars);

    // gather set of ids which are hit
    for (unsigned i = 0; i < nObjVars; ++i) {
      if (IloAbs(vals[i] - 1.0) < EPS) {
        int id;
        sscanf(objVars[i].getName(), "%10d", &id);
        out_hittingSet.insert(id);
      }
    }

    // TODO: recompute integer weight
    LB = std::max(LB, (unsigned long)cplex.getObjValue());

  } catch (IloException e) {
    cerr << e << endl;

    exit(1);
  }

  return true;
}

uint64_t CPLEXInterface::computeRelaxedBound(
    uint64_t & LB, uint64_t UB, unordered_map<int,bool> &forcedVars)
{

  double LB_lp = 0;

  try {

    // attach objective function if detached
    if (!objFuncAttached) {
      model.add(objective);
      objFuncAttached = true;
    }

    if (lp_cplex.isMIP()) {
      log << "c adding LP conversions" << endl;
      lp_model.add(IloConversion(*env, vars, ILOFLOAT));
      lp_model.add(IloConversion(*env, objVars, ILOFLOAT));
    }

    if (lp_cplex.isMIP()) {
      cerr << "c LP relaxation failed" << endl;
      exit(1);
    }

    IloNumArray reducedCosts(*env);
    IloNumArray lp_vals(*env);

    try {
      ++lpCalls;
      lpTimer.start();
      lp_cplex.solve();
      lpTimer.stop();

      LB_lp = uint64_t(lp_cplex.getObjValue());
      log << "c LP LB " << LB_lp << endl;

      if (UB == int64_t(ceil(LB_lp))) {
        LB = UB;
        return LB;
      }

      lp_cplex.getReducedCosts(reducedCosts, objVars);
      lp_cplex.getValues(lp_vals, objVars);

    } catch (IloCplex::Exception e) {
      cerr << e.getMessage() << endl;
      exit(1);
    }

    //
    // Reduced cost fixing
    //
    if (cfg.reducedCosts)
    for (unsigned i = 0; i < nObjVars; ++i) {
      // variable is forced already?
      if (objVars[i].getUB() == objVars[i].getLB())
        continue;

      int bVar;
      sscanf(objVars[i].getName(), "%10d", &bVar);
      int64_t w = var_to_weight[bVar];

      double rc = reducedCosts[i];
      int64_t bvar_LB = 0;

      // conditional lb cannot exceed ub
      if (UB >= int64_t(ceil(LB_lp + w - EPS))) {
        continue;
      }

      //assert(instance->UB_bool_solution.size() > unsigned(bVar));

      if (IloAbs(lp_vals[i]) < EPS) { // try to harden

        // reduced cost is zero or negative?
        if (rc <= EPS) continue;

        bvar_LB = int64_t(ceil(LB_lp + rc - EPS));

        //log << "c bvar_LB=" << bvar_LB << "\t" << " UB=" << UB << endl;

        // conditinal lb too low to force
        if (UB > bvar_LB) continue;

        //bool true_in_best_model = instance->UB_bool_solution[bVar];
        // cannot force to 0 since best found model might be optimal
        if (/*true_in_best_model && */bvar_LB == UB) continue;

        //forceVar(instance, bVar, false);
        objVars[i].setBounds(0, 0);
        forcedVars[bVar] = 0;

        log << "c forced bVar " << bVar << " lp_lb: " << LB_lp \
          << " ub: " << UB << " bv_lb: " << bvar_LB << endl;
        ++reducedCostForcedVars;

      } else {  // try to relax

        // reduced cost is zero or positive?
        if (rc >= -EPS) continue;

        bvar_LB = int64_t(ceil(LB_lp - rc - EPS));

        //log << "c bvar_LB=" << bvar_LB << "\t" << " UB=" << UB << endl;

        // conditinal lb too low to relax
        if (UB > bvar_LB) continue;

        //bool false_in_best_model = !instance->UB_bool_solution[bVar];
        // cannot force to 1 since best found model might be optimal
        if (/*false_in_best_model && */bvar_LB == UB) continue;

        //forceVar(instance, bVar, true);
        objVars[i].setBounds(1, 1);
        forcedVars[bVar] = 1;

        log << "c relaxed bVar " << bVar << " lp_lb: " << LB_lp \
          << " ub: " << UB << " bv_lb: " << bvar_LB << endl;
        ++reducedCostRelaxedVars;
      }
    }

  } catch (IloException e) {
    cerr << e << endl;
    exit(1);
  }

  return int64_t(ceil(LB_lp - EPS));
}

void CPLEXInterface::printStats() {
  log << "c IP calls " << ipCalls << endl;
  log << "c IP time " << ipTimer.cpu_ms_total() << " ms" << endl;
  if (cfg.reducedCosts) {
    log << "c LP calls " << lpCalls << endl;
    log << "c LP time " << lpTimer.cpu_ms_total() << endl;
    log << "c rc fixes " << reducedCostRelaxedVars + reducedCostForcedVars << endl;
  }
}
