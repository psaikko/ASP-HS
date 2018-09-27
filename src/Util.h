#pragma once

#include <string>
#include <vector>
#include <limits.h> 		 // INT_MAX etc
#include <cstdarg>
#include <ostream>

template <typename T>
inline std::ostream& operator<<(std::ostream& o, const std::vector<T> & v) {
	if (v.size())
		o << v[0];
	for (unsigned i = 1; i < v.size(); ++i)
		o << " " << v[i];
	return o;
}

bool string_contains(const std::string & a, const std::string & b);

bool string_beginsWith(const std::string & a, const std::string & b);

bool string_endsWith(const std::string & a, const std::string & b);

// a subset of b ?
bool core_subset(const std::vector<int> & a, const std::vector<int> & b);

std::string trim(const std::string & s);

void terminate(const int code, const char * fmt, ...);

void condTerminate(const bool cond, const int code, const char * fmt, ...);

void log(const int level, const char * fmt, ...);

void condLog(const bool cond, const int level, const char * fmt, ...);

void logCore(const int level, std::vector<int> & core);

#ifndef ASPHS_DEBUG

#define traceMsg
#define traceIf
#define traceAction

#else

#define traceMsg(level, message) \
do { \
  if (level <= GlobalConfig::get().traceLevel) { \
    for (int _log_indent = 0; _log_indent < level; ++_log_indent) std::cerr << "  "; \
    std::cerr << message << std::endl; \
  } \
} while (0)

#define traceIf(level, condition, message) \
do { \
  if (condition) { \
    traceMsg(level, message) \
  } \
} while (0)

#define traceAction(level, action) \
do { \
  if (level <= GlobalConfig::get().traceLevel) { \
    action \
  } \
} while (0)

#endif