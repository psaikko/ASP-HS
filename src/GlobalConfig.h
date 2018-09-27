#pragma once

#include <vector>
#include <utility>       // pair
#include <unordered_map>
#include <string>
#include <iosfwd>

class GlobalConfig {
 public:

  static inline GlobalConfig& get() {
    static GlobalConfig config;
    return config;
  }

  bool initialized;
  bool stdin;

  #include "OptionDeclarations.cpp"

  void printHelp(std::ostream & out);
  void parseArgs(int argc, const char** argv, std::ostream & out);

 private:
  GlobalConfig() { initialized = 0; stdin = 0; };
  GlobalConfig(GlobalConfig const&);
  void operator=(GlobalConfig const&);

};
