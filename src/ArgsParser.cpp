#include "ArgsParser.h"
#include "Util.h"

#include <algorithm>

using namespace std;

bool ArgsParser::hasFlag(const string& option) {
  return find(argsBegin, argsEnd, option) != argsEnd;
}

bool ArgsParser::getBoolOption(const string& option, bool defaultValue) {
  if (hasFlag("--no-" + option)) return false;
  if (hasFlag("--" + option))
    return true;
  else
    return defaultValue;
}

long int ArgsParser::getIntOption(const string& option, long int defaultValue) {
  if (!hasFlag("--" + option)) return defaultValue;
  long int value = defaultValue;
  if (!sscanf(getOption("--" + option), "%100ld", &value))
    terminate(1, "--%s expects parameter\n", option.c_str());
  return value;
}

double ArgsParser::getDoubleOption(const string& option, double defaultValue) {
  if (!hasFlag("--" + option)) return defaultValue;
  double value = defaultValue;
  if (!sscanf(getOption("--" + option), "%100lf", &value))
    terminate(1, "--%s expects parameter\n", option.c_str());
  return value;
}

string ArgsParser::getStringOption(const string& option,
                                   const string& defaultValue) {
  if (!hasFlag("--" + option)) return defaultValue;
  const char* value = getOption("--" + option);
  if (!value || (value[0] == '-' && value[1] == '-')) 
    terminate(1, "--%s expects parameter\n", option.c_str());
  return trim(string(value));
}

const char* ArgsParser::getOption(const string& option) {
  const char** itr = find(argsBegin, argsEnd, option);
  if (itr != argsEnd && ++itr != argsEnd) return *itr;
  return "";
}