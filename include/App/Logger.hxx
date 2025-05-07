#ifndef LOGGER_HXX
#define LOGGER_HXX

#include <iostream>
#include <sstream>
#include <string>

#define INFO(fmt, ...)                            \
    do {                                          \
        std::ostringstream oss;                   \
        oss << "[INFO ] " << fmt;                 \
        printf(oss.str().c_str(), ##__VA_ARGS__); \
        std::cout << '\n';                        \
    } while (0)

#define WARNING(fmt, ...)                         \
    do {                                          \
        std::ostringstream oss;                   \
        oss << "[WARN ] " << fmt;                 \
        printf(oss.str().c_str(), ##__VA_ARGS__); \
        std::cout << '\n';                        \
    } while (0)

#define ERROR(fmt, ...)                                    \
    do {                                                   \
        std::ostringstream oss;                            \
        oss << "[ERROR] " << fmt;                          \
        fprintf(stderr, oss.str().c_str(), ##__VA_ARGS__); \
        std::cerr << '\n';                                 \
    } while (0)

#endif  // LOGGER_HXX
