#ifndef TIMER_H
#define TIMER_H

#include "Profiler.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>

class Timer
{
public:
    Timer(std::string name) : _name(name)
    {
        _start_time = std::chrono::high_resolution_clock::now();
    }

    ~Timer()
    {
        double elapsed_time = (std::chrono::high_resolution_clock::now() - _start_time).count();
        Profiler::instance().add(_name, elapsed_time);
    }

private:
    std::string _name;
    std::chrono::high_resolution_clock::time_point _start_time;
};

#endif // TIMER_H