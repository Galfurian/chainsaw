/// @file it_algebra.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Provides efficient algebraic operations such as computing sums,
/// differences, and scaled sums over iterators. Simplification of the code
/// available at:
///     https://github.com/headmyshoulder/odeint-v2
/// I've also added some more functions.

#pragma once

#include <cmath>
#include <limits>
#include <functional>

namespace numint::detail::it_algebra
{

/// @brief Computes the maximum absolute difference between elements in two ranges.
/// @param a0_first Iterator to the first element of range 1.
/// @param a0_last Iterator to the last element of range 1.
/// @param a1_first Iterator to the first element of range 2.
/// @param a1_last Iterator to the last element of range 2.
/// @return Maximum absolute difference.
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

/// @brief Computes the maximum relative difference between elements in two ranges.
/// @param a0_first Iterator to the first element of range 1.
/// @param a0_last Iterator to the last element of range 1.
/// @param a1_first Iterator to the first element of range 2.
/// @param a1_last Iterator to the last element of range 2.
/// @return Maximum relative difference.
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

/// @brief Computes the maximum of absolute and relative differences between two ranges.
/// @param a0_first Iterator to the first element of range 1.
/// @param a0_last Iterator to the last element of range 1.
/// @param a1_first Iterator to the first element of range 2.
/// @param a1_last Iterator to the last element of range 2.
/// @return Maximum combined absolute and relative difference.
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

namespace detail
{

/// @brief Base case for the recursive variadic function that adds scaled terms.
/// @param ... Unused parameters for recursion termination.
/// @note When there are no more scalars or iterators, the recursion stops.
constexpr inline void add_helper(...) noexcept
{
    // Base case: Do nothing, recursion stops here.
}

/// @brief Helper function to recursively add scaled terms without modifying iterators prematurely.
/// @param y The value to accumulate the scaled terms into.
/// @param op Operation to perform on the scalar and dereferenced iterator value.
/// @param a Scalar value to scale the current term.
/// @param x Iterator corresponding to the scalar.
/// @param args Remaining scalars and iterators.
/// @note This function is called recursively to handle multiple terms.
template <class T, class It, class Op, class... Args>
constexpr inline void add_helper(T &y, Op op, T a, It &x, Args &...args) noexcept
{
    // Add the current scaled term.
    y += op(a, *x++);
    // Recursively process the remaining terms.
    add_helper(y, op, args...);
}

} // namespace detail

/// @brief Computes the element-wise sum of multiple scaled ranges into the output range.
/// @param y_first Iterator to the first element of the output range.
/// @param y_last Iterator to the last element of the output range.
/// @param op Operation to apply for element-wise computation (e.g., addition, subtraction).
/// @param a First scalar value to scale the first range.
/// @param x Iterator corresponding to the first scalar.
/// @param args Variadic template to accept additional scalars and iterators for further ranges.
/// @note This function uses variadic templates to accept any number of scalars and corresponding iterators.
template <class OutIt, class T, class InIt, class Op, class... Args>
constexpr inline void sum_operation(OutIt y_first, OutIt y_last, Op op, T a, InIt x, Args... args) noexcept
{
    while (y_first != y_last) {
        // Add the current scaled term.
        *y_first = op(a, *x++);
        // Recursively add the remaining scalars and iterators.
        detail::add_helper(*y_first, op, args...);
        // Increment the output iterator.
        ++y_first;
    }
}

/// @brief Accumulates the element-wise sum of multiple scaled ranges into the output range.
/// @param y_first Iterator to the first element of the output range.
/// @param y_last Iterator to the last element of the output range.
/// @param op Operation to apply for element-wise computation (e.g., addition, subtraction).
/// @param a First scalar value to scale the first range.
/// @param x Iterator corresponding to the first scalar.
/// @param args Variadic template to accept additional scalars and iterators for further ranges.
/// @note This function uses variadic templates to accept any number of scalars and corresponding iterators.
template <class OutIt, class T, class InIt, class Op, class... Args>
constexpr inline void accumulate_operation(OutIt y_first, OutIt y_last, Op op, T a, InIt x, Args... args) noexcept
{
    while (y_first != y_last) {
        // Add the current scaled term.
        *y_first += op(a, *x++);
        // Recursively add the remaining scalars and iterators.
        detail::add_helper(*y_first, op, args...);
        // Increment the output iterator.
        ++y_first;
    }
}

} // namespace numint::detail::it_algebra
