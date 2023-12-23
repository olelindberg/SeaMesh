#include "Profiler.h"
#include <iomanip>

bool timings_cmp(std::pair<std::string, double> &a, std::pair<std::string, double> &b)
{
    return a.second > b.second;
}

Profiler &Profiler::instance()
{
    static Profiler instance;
    return instance;
}

void Profiler::add(std::string name, double time)
{
    if (_timings.find(name) == _timings.end())
        _timings.insert(std::make_pair(name, time));
    else
        _timings[name] += time;
}

void Profiler::print()
{
    std::cout << "Printing timing results ..." << std::endl;

    std::vector<std::pair<std::string, double>> sorted_timings;
    for (auto &it : _timings)
        sorted_timings.push_back(it);
    std::sort(sorted_timings.begin(), sorted_timings.end(), timings_cmp);
    for (auto &it : sorted_timings)
    {

        std::cout << std::left << std::setw(25) << it.first << std::right << " " << it.second << std::endl;
    }

    std::cout << "Printing timing results, done" << std::endl;
}