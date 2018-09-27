#include "Timer.h"
#include <cassert>

using namespace std::chrono;

void Timer::start() {
	start_clocks = clock();
	start_time  = system_clock::now();
	assert(!counting);
	counting = true;
}

void Timer::stop() {
	clock_t end_clocks = clock();
	auto end_time = system_clock::now();
	assert(counting);
	counting = false;

	total_clocks += (end_clocks - start_clocks);
	total_us += duration_cast<microseconds>(end_time - start_time).count();
}

unsigned long Timer::real_ms_current() {
	if (!counting) return 0;
	auto current_time = system_clock::now();
	return duration_cast<microseconds>(current_time - start_time).count() / 1000;
}

unsigned long Timer::cpu_ms_current() {
	if (!counting) return 0;
	clock_t current_clocks = clock();
	return 1000 * (current_clocks - start_clocks) / CLOCKS_PER_SEC;
}

unsigned long Timer::real_ms_total() {
	return total_us / 1000 + real_ms_current();
}

unsigned long Timer::cpu_ms_total() {
	return 1000 * total_clocks / CLOCKS_PER_SEC + cpu_ms_current();
}

void Timer::add(Timer &other) {
	assert(!other.counting);

	total_clocks += other.total_clocks;
	total_us += other.total_us;
}