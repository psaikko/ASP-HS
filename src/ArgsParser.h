#pragma once

#include <string>

class ArgsParser {
 public:
  ArgsParser(const char** begin, const char** end) : argsBegin(begin), argsEnd(end) { }

  bool hasFlag(const std::string& option);
  bool getBoolOption(const std::string& option, bool defaultValue);
  long int getIntOption(const std::string& option, long int defaultValue);
  double getDoubleOption(const std::string& option, double defaultValue);
  std::string getStringOption(const std::string& option,
                              const std::string& defaultValue);

 protected:
  const char* getOption(const std::string& option);
  const char** argsBegin, **argsEnd;
};
