#pragma once

#include <ilcplex/ilocplex.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iosfwd>
#include <cmath>

#include "Timer.h"
#include "GlobalConfig.h"

class CPLEXInterface {

 public:
  CPLEXInterface(std::ostream & log = std::cout);
  ~CPLEXInterface();

  void addObjectiveVariable(int bVar, uint64_t weight);

  void addVariable(int var);

  void addCore(std::vector<int>& core);

  void addClause(std::vector<int>& cl, bool lp_only = false);

  void addAggregate(int id, std::vector<int> literals, std::vector<uint64_t> weights, uint64_t bound, bool lp_only = false);

  bool computeHS(std::unordered_set<int>& out_hittingSet,
    uint64_t & LB, uint64_t UB, std::unordered_map<int,bool> &forcedVars);

  uint64_t computeRelaxedBound(
    uint64_t & LB, uint64_t UB, std::unordered_map<int,bool> &forcedVars);

  void printStats();

 private:

  std::ostream & log;
  GlobalConfig & cfg;

  IloNumVar newObjVar(int v);

  bool objFuncAttached;

  IloEnv* env;
  IloModel model;
  IloModel lp_model;
  IloNumVarArray objVars;
  IloNumVarArray vars;
  unsigned nObjVars;
  unsigned nVars;
  IloObjective objective;
  IloRangeArray cons;
  IloCplex cplex;
  IloCplex lp_cplex;
  std::unordered_map<int, IloNumVar> var_to_IloVar;
  std::unordered_map<int, uint64_t> var_to_weight;

  Timer ipTimer;
  Timer lpTimer;

  unsigned ipCalls;
  unsigned lpCalls;

  unsigned reducedCostForcedVars;
  unsigned reducedCostRelaxedVars;
};
