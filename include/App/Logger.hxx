#ifndef T2S_LOGGER_HXX
#define T2S_LOGGER_HXX

#include <cstdio>
#include <utility>

#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/format.h>

namespace Logger {

template <typename... Args>
inline void Debug(std::string_view func, fmt::format_string<Args...> format, Args&&... args) {
    fmt::print(stdout, "{}", fmt::styled("DEBUG", fmt::fg(fmt::color::white) | fmt::emphasis::bold));
    fmt::print(stdout, " :: {} :: ", func);
    fmt::print(stdout, format, std::forward<Args>(args)...);
    fmt::print(stdout, "\n");
}

template <typename... Args>
inline void Info(std::string_view func, fmt::format_string<Args...> format, Args&&... args) {
    fmt::print(stdout, "{}", fmt::styled("INFO ", fmt::fg(fmt::color::green) | fmt::emphasis::bold));
    fmt::print(stdout, " :: {} :: ", func);
    fmt::print(stdout, format, std::forward<Args>(args)...);
    fmt::print(stdout, "\n");
}

template <typename... Args>
inline void Warning(std::string_view func, fmt::format_string<Args...> format, Args&&... args) {
    fmt::print(stderr, "{}", fmt::styled("WARN ", fmt::fg(fmt::color::yellow) | fmt::emphasis::bold));
    fmt::print(stderr, " :: {} :: ", func);
    fmt::print(stderr, format, std::forward<Args>(args)...);
    fmt::print(stderr, "\n");
}

template <typename... Args>
inline void Error(std::string_view func, fmt::format_string<Args...> format, Args&&... args) {
    fmt::print(stderr, "{}", fmt::styled("ERROR", fmt::fg(fmt::color::red) | fmt::emphasis::bold));
    fmt::print(stderr, " :: {} :: ", func);
    fmt::print(stderr, format, std::forward<Args>(args)...);
    fmt::print(stderr, "\n");
}

}  // namespace Logger

#endif  // T2S_LOGGER_HXX
