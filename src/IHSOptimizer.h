#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <cmath>

#include "Timer.h"
#include "GlobalConfig.h"

#include "ProgramListener.h"


class CPLEXInterface;
class WaspFacade;

struct OptimizationLiteralData;
struct Literal;

class IHSOptimizer : public ProgramListener {
 public:
  IHSOptimizer(std::ostream & log = std::cout);
  ~IHSOptimizer();

  GlobalConfig & cfg;

  std::ostream & log;

  std::vector< OptimizationLiteralData > optLits;

  std::vector< std::vector< int > > cores;

  std::unordered_map< int, uint64_t > id_weight;
  std::unordered_map< int, uint64_t > id_coreCount;


  void solve(std::istream & instance);

  void setAssumptions(std::unordered_set< int > & hittingSet, std::vector< Literal > & assumptions);

  uint64_t answerSetCost();

  void foundCore(std::vector< Literal > & core);

  void updateUB(uint64_t cost);

  void minimizeCore(std::vector< Literal > & core);

  void cardMinimize(std::vector< Literal > & core);

  Literal addTemporaryCard(std::vector< Literal > & core, unsigned int bound);

  void eqSeeding();

  void feedClauses(CPLEXInterface * solver, bool lp_only = false);

  std::unordered_set< int > greedyHittingSet();

  void printStats();

  uint64_t UB;
  uint64_t LB;

  Timer solveTimer;
  Timer coreExtractionTimer;
  Timer minimizeTimer;

  double avgCoreSize;
  double avgCoreSizeMinimized;

  void addedClause(const std::vector<int>& literals);

  void addedAggregate(unsigned int id,
    const std::vector<int>& literals,
    const std::vector<uint64_t>& weights,
    uint64_t bound);

 private:
  WaspFacade * ASP_solver;
  CPLEXInterface * IP_solver;

  std::vector<std::vector<int>> clauses;

  typedef struct Aggregate {
    unsigned int id;
    std::vector<int> literals;
    std::vector<uint64_t> weights;
    uint64_t bound;
  } Aggregate;

  std::vector<Aggregate> aggregates;
};
