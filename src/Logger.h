#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <ctime>
#include <iostream>
#include <mutex>
#include <string>

class Logger {
public:
  // Log an info message
  static void info(const std::string &msg) { log("INFO", msg); }

  // Log an error message
  static void error(const std::string &msg) { log("ERROR", msg); }

private:
  // Thread-safe logging
  static void log(const std::string &level, const std::string &msg) {
    static std::mutex           mtx;
    std::lock_guard<std::mutex> lock(mtx);

    // Timestamp
    auto        now = std::chrono::system_clock::now();
    std::time_t t   = std::chrono::system_clock::to_time_t(now);

    // Trim newline from ctime
    std::string ts = std::ctime(&t);
    ts.pop_back();

    std::cout << "[" << ts << "] [" << level << "] " << msg << std::endl;
  }
};

#endif // LOGGER_H