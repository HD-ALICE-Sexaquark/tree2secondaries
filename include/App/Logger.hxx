#pragma once

#include <cstdio>
#include <format>
#include <print>
#include <string_view>
#include <utility>

namespace Logger {

template <typename... Args>
inline void Debug(std::string_view func, std::format_string<Args...> fmt, Args &&...args) {
    std::println(stdout, "DEBUG :: {} :: {}", func, std::format(fmt, std::forward<Args>(args)...));
}

template <typename... Args>
inline void Info(std::string_view func, std::format_string<Args...> fmt, Args &&...args) {
    std::println(stdout, "INFO  :: {} :: {}", func, std::format(fmt, std::forward<Args>(args)...));
}

template <typename... Args>
inline void Warning(std::string_view func, std::format_string<Args...> fmt, Args &&...args) {
    std::println(stderr, "WARN  :: {} :: {}", func, std::format(fmt, std::forward<Args>(args)...));
}

template <typename... Args>
inline void Error(std::string_view func, std::format_string<Args...> fmt, Args &&...args) {
    std::println(stderr, "ERROR :: {} :: {}", func, std::format(fmt, std::forward<Args>(args)...));
}

}  // namespace Logger
