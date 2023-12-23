#ifndef PROFILER_H
#define PROFILER_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>

class Profiler
{
public:
    Profiler(Profiler const &) = delete;
    void operator=(Profiler const &) = delete;

    static Profiler &instance();

    void add(std::string name, double time);

    void print();

private:
    std::map<std::string, double> _timings;

    Profiler() {}
};

#endif // PROFILER_H