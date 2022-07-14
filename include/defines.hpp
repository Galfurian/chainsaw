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

template <typename T>
constexpr inline T compute_samples(Time time_start, Time time_end, Time time_delta, Time sampling = 1.0)
{
    return static_cast<T>(((time_end - time_start) / time_delta) * sampling);
}

/// @brief Returns the sign of the input value (+1 or -1).
template <typename T>
constexpr inline int sign(const T &value)
{
    return (value > 0) - (value < 0);
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(std::array<T1, N> a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] += b;
    return a64l;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator-(std::array<T1, N> a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] -= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*(std::array<T1, N> a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/(std::array<T1, N> a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] /= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(T1 a, std::array<T2, N> b)
{
    for (std::size_t i = 0; i < N; ++i)
        b[i] += a;
    return b;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator-(T1 a, std::array<T2, N> b)
{
    for (std::size_t i = 0; i < N; ++i)
        b[i] -= a;
    return b;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*(T1 a, std::array<T2, N> b)
{
    for (std::size_t i = 0; i < N; ++i)
        b[i] *= a;
    return b;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/(T1 a, std::array<T2, N> b)
{
    for (std::size_t i = 0; i < N; ++i)
        b[i] /= a;
    return b;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+=(std::array<T1, N> &a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] += b;
    return a64l;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator-=(std::array<T1, N> &a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] -= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*=(std::array<T1, N> &a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/=(std::array<T1, N> &a, T2 b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] /= b;
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+(std::array<T1, N> a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] += b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator-(std::array<T1, N> a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] -= b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*(std::array<T1, N> a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/(std::array<T1, N> a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] /= b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator+=(std::array<T1, N> &a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] += b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator-=(std::array<T1, N> &a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] -= b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator*=(std::array<T1, N> &a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b[i];
    return a;
}

template <class T1, class T2, std::size_t N>
constexpr inline auto operator/=(std::array<T1, N> &a, const std::array<T2, N> &b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] /= b[i];
    return a;
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
