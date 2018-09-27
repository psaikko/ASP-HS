#pragma once

#include <ctime>
#include <chrono>

class Timer {
 public:
	Timer() : total_us(0), total_clocks(0), counting(false) {};

	void start();
	void stop();

	unsigned long real_ms_total();
	unsigned long real_ms_current();

	unsigned long cpu_ms_total();
	unsigned long cpu_ms_current();

	void add(Timer &other);

 protected:

 	std::chrono::time_point<std::chrono::system_clock> start_time;
 	unsigned long total_us;

 	clock_t total_clocks;
 	clock_t start_clocks;

 	bool counting;
};