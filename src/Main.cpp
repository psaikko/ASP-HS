#include "IHSOptimizer.h"
#include "GlobalConfig.h"

#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include "util/WaspConstants.h"
#include "util/WaspOptions.h"

using namespace std;

int EXIT_CODE = 0;

IHSOptimizer * solver;

void stop(int) {
  if (solver)
    solver->printStats();
  printf("s UNKNOWN\n");
  fflush(stdout);
  _Exit(1);
}

int main(int argc, const char* argv[]) {

  GlobalConfig & cfg = GlobalConfig::get();

  clock_t begin = clock();

  if(argc < 2) {
    printf("Error: need filepath.\n");
    exit(1);
  }

  cfg.parseArgs(argc, argv, cout);

  signal(SIGINT, stop);
  signal(SIGTERM, stop);
  signal(SIGXCPU, stop);

  ofstream nullstream("/dev/null");
  ostream & log = cfg.silent ? nullstream : cout;

  solver = new IHSOptimizer(log);

  if (cfg.stdin) {
    solver->solve(cin);
  } else {
    ifstream file(argv[1]);
    solver->solve(file);
  }

  solver->printStats();
  log << "c CPU time: " << (float(clock() - begin) / float(CLOCKS_PER_SEC)) << " seconds" << endl;
  return 0;
}
