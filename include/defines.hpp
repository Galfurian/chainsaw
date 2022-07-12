/// @file defines.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-04-13

#pragma once

#include <iostream>
#include <array>

/// The time is a continuous time value.
using Time = double;

/// State variables are declared as floating point.
using Variable = double;

/// @brief Returns the sign of the input value (+1 or -1).
template <typename T>
constexpr inline int sign(const T &value)
{
    return (value > 0) - (value < 0);
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(const std::array<T1, N> &a, T2 b)
{
    std::array<T1, N> result(a);
    for (auto &i : result)
        i += b;
    return result;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(T1 a, const std::array<T2, N> &b)
{
    return b + a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(const std::array<T1, N> &a, const std::array<T2, N> &b)
{
    std::array<T1, N> result(a);
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] += b[i];
    return result;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+=(std::array<T1, N> &a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] += b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*(const std::array<T1, N> &a, T2 b)
{
    std::array<T1, N> result(a);
    for (auto &i : result)
        i *= b;
    return result;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*(T1 a, const std::array<T2, N> &b)
{
    return b * a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/(const std::array<T1, N> &a, T2 b)
{
    std::array<T1, N> result(a);
    for (auto &i : result)
        i /= b;
    return result;
}

template <class T, std::size_t N>
T abs(const std::array<T, N> &a)
{
    T result = 0;
    for (const auto &i : a)
        result += std::abs(i);
    return result;
}

template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &lhs, const std::array<T, N> &rhs)
{
    lhs << "[ ";
    for (std::size_t i = 0; i < N; ++i)
        std::cout << rhs[i] << " ";
    lhs << "]";
    return lhs;
}
