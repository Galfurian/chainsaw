/// @file it_algebra.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2
/// I've also added some more functions.

#pragma once

#include <cmath>
#include <limits>

namespace solver::detail::it_algebra
{

/// @brief Computes the maximum absolute difference as : sum(abs(a1 - a2))
/// @param a0_first
/// @param a0_last
/// @param a1_first
/// @param a1_last
/// @return constexpr T
template <class T, class It>
constexpr inline T max_abs_diff(It a0_first, It a0_last, It a1_first, It a1_last) noexcept
{
    // Initialize the value to epsilon, to prevent small truncation error when
    // using the returned value.
    T ret(std::numeric_limits<T>::epsilon());
    // Find the max difference.
    while ((a0_first != a0_last) && (a1_first != a1_last)) {
        ret = std::max(ret, std::abs(*a0_first - *a1_first));
        ++a0_first, ++a1_first;
    }
    return ret;
}

/// @brief Computes the maximum relative difference as : sum(abs(a1 - a2) / a1)
/// @param a0_first
/// @param a0_last
/// @param a1_first
/// @param a1_last
/// @return constexpr T
template <class T, class It>
constexpr inline T max_rel_diff(It a0_first, It a0_last, It a1_first, It a1_last) noexcept
{
    // Initialize the value to epsilon, to prevent small truncation error when
    // using the returned value.
    T ret(std::numeric_limits<T>::epsilon());
    // Find the max difference.
    while ((a0_first != a0_last) && (a1_first != a1_last)) {
        ret = std::max(ret, std::abs((*a0_first - *a1_first) / *a0_first));
        ++a0_first, ++a1_first;
    }
    return ret;
}

/// @brief Computes the maximum between absolute and relative difference.
/// @param a0_first
/// @param a0_last
/// @param a1_first
/// @param a1_last
/// @return constexpr T
template <class T, class It>
constexpr inline T max_comb_diff(It a0_first, It a0_last, It a1_first, It a1_last) noexcept
{
    // Initialize the value to epsilon, to prevent small truncation error when
    // using the returned value.
    T ret(std::numeric_limits<T>::epsilon());
    // Find the max difference.
    while ((a0_first != a0_last) && (a1_first != a1_last)) {
        ret = std::max(std::max(ret, std::abs((*a0_first - *a1_first) / *a0_first)), std::abs(*a0_first - *a1_first));
        ++a0_first, ++a1_first;
    }
    return ret;
}

// computes sum(y)
template <class T, class It>
constexpr inline T accumulate(It y_first, It y_last) noexcept
{
    T ret = T(0);
    while (y_first != y_last) {
        ret += (*y_first++);
    }
    return ret;
}

// computes y = abs(y)
template <class It>
constexpr inline void abs(It y_first, It y_last) noexcept
{
    while (y_first != y_last) {
        (*y_first) = std::abs(*y_first);
        ++y_first;
    }
}

// computes sum(abs(y))
template <class T, class It>
constexpr inline T accumulate_abs(It y_first, It y_last) noexcept
{
    T ret = T(0);
    while (y_first != y_last) {
        ret += std::abs(*y_first++);
    }
    return ret;
}

// computes y = a1 * x1
template <class OutIt, class InIt, class T>
constexpr inline void increment(OutIt y_first, OutIt y_last, InIt x1, T a) noexcept
{
    while (y_first != y_last) {
        (*y_first++) += a * (*x1++);
    }
}

// computes y = x1 + x2
template <class OutIt, class InIt1, class InIt2>
constexpr inline void sum(OutIt y_first, OutIt y_last, InIt1 x1, InIt2 x2) noexcept
{
    while (y_first != y_last) {
        (*y_first++) = (*x1++) + (*x2++);
    }
}

// computes y = x1 - x2
template <class OutIt, class InIt1, class InIt2>
constexpr inline void sub(OutIt y_first, OutIt y_last, InIt1 x1, InIt2 x2) noexcept
{
    while (y_first != y_last) {
        (*y_first++) = (*x1++) - (*x2++);
    }
}

// computes y = a1*x1 + a2*x2
template <class OutIt, class InIt1, class InIt2, class T>
constexpr inline void scale_sum(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2) noexcept
{
    while (y_first != y_last) {
        (*y_first++) = a1 * (*x1++) + a2 * (*x2++);
    }
}

// computes y = x1 + a*x2
template <class OutIt, class InIt1, class InIt2, class T>
constexpr inline void scale_sum(OutIt y_first, OutIt y_last, InIt1 x1, T a, InIt2 x2) noexcept
{
    while (y_first != y_last) {
        (*y_first++) = (*x1++) + a * (*x2++);
    }
}

// computes y = x1 + a2*x2 + a3*x3
template <class OutIt, class InIt1, class InIt2, class InIt3, class T>
constexpr inline void scale_sum(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2, T a3, InIt3 x3) noexcept
{
    while (y_first != y_last) {
        (*y_first++) = a1 * (*x1++) + a2 * (*x2++) + a3 * (*x3++);
    }
}

// computes y += a1*x1 + a2*x2 + a3*x3 + a4*x4
template <class OutIt, class InIt1, class InIt2, class InIt3, class InIt4, class T>
constexpr inline void scale_sum_inplace(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2, T a3, InIt3 x3, T a4, InIt4 x4) noexcept
{
    while (y_first != y_last) {
        (*y_first++) += a1 * (*x1++) + a2 * (*x2++) + a3 * (*x3++) + a4 * (*x4++);
    }
}

// computes tmp = y, y = x1 + a*x2, x1 = tmp
template <class OutIt, class InIt, class T>
constexpr inline void scale_sum_swap(OutIt y_first, OutIt y_last, OutIt x1, T a, InIt x2) noexcept
{
    T swap = static_cast<T>(.0);
    while (y_first != y_last) {
        swap       = (*x1) + a * (*x2++);
        *x1++      = *y_first;
        *y_first++ = swap;
    }
}

} // namespace solver::detail::it_algebra
