#include <string>
#include <ostream>
#include <climits>

#include "GlobalConfig.h"
#include "ArgsParser.h"
#include "Util.h"

#include "OptionParsing.h"

using namespace std;

void GlobalConfig::printHelp(ostream & out) {
  
  #include "OptionHelp.cpp"

  out << "For integer/double/string -valued options use \"--option value\"" << endl;
  out << "For boolean-valued options use \"--option\" or \"--no-option\"" << endl;

  exit(0);
}

// handle user-settable options from command line arguments
void GlobalConfig::parseArgs(int argc, const char** argv, ostream & out) {

  ArgsParser args(argv, argv + argc);

  unordered_map<string, string> flag_type;

  #include "OptionParsing.cpp"

  // validate flags
  if (argv[1] == string("-")) {
    stdin = true;
  }

  for (int i = 2; i < argc; ++i) {
    for (auto p : flag_type) {
      string flag = p.first;
      string type = p.second;
      
      if (argv[i] == ("--" + flag)) {
        if (type == "bool") {
          goto next;
        } else {
          // skip parameter
          goto skip;
        }
      }

      if (type == "bool") {
        if (argv[i] == ("--no-" + flag)) {
          goto next;
        }
      }
    }

    cout << "Unexpected parameter " << argv[i] << endl;
    exit(1);

    skip: i += 1;

    next: continue;
  }

  initialized = true;

  if (help != 0) {
    printHelp(out);
  }
}
