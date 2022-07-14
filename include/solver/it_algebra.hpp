/// @file it_algebra.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

namespace solver::it_algebra
{

// computes sum(y)
template <class T, class OutIt>
constexpr inline T accumulate(OutIt y_first, OutIt y_last) noexcept
{
    T ret = T(0);
    while (y_first != y_last)
        ret += (*y_first++);
    return ret;
}

// computes y = abs(y)
template <class OutIt>
constexpr inline void abs(OutIt y_first, OutIt y_last) noexcept
{
    while (y_first != y_last) {
        (*y_first) = std::abs(*y_first);
        ++y_first;
    }
}

// computes sum(abs(y))
template <class T, class OutIt>
constexpr inline T accumulate_abs(OutIt y_first, OutIt y_last) noexcept
{
    T ret = T(0);
    while (y_first != y_last)
        ret += std::abs(*y_first++);
    return ret;
}

// computes y = a1 * x1
template <class OutIt, class InIt, class T>
constexpr inline void increment(OutIt y_first, OutIt y_last, InIt x1, T a) noexcept
{
    while (y_first != y_last)
        (*y_first++) += a * (*x1++);
}

// computes y = x1 + x2
template <class OutIt, class InIt1, class InIt2>
constexpr inline void sum(OutIt y_first, OutIt y_last, InIt1 x1, InIt2 x2) noexcept
{
    while (y_first != y_last)
        (*y_first++) = (*x1++) + (*x2++);
}

// computes y = x1 - x2
template <class OutIt, class InIt1, class InIt2>
constexpr inline void sub(OutIt y_first, OutIt y_last, InIt1 x1, InIt2 x2) noexcept
{
    while (y_first != y_last)
        (*y_first++) = (*x1++) - (*x2++);
}

// computes y = a1*x1 + a2*x2
template <class OutIt, class InIt1, class InIt2, class T>
constexpr inline void scale_sum(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2) noexcept
{
    while (y_first != y_last)
        (*y_first++) = a1 * (*x1++) + a2 * (*x2++);
}

// computes y = x1 + a2*x2 + a3*x3
template <class OutIt, class InIt1, class InIt2, class InIt3, class T>
constexpr inline void scale_sum(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2, T a3, InIt3 x3) noexcept
{
    while (y_first != y_last)
        (*y_first++) = a1 * (*x1++) + a2 * (*x2++) + a3 * (*x3++);
}

// computes y += a1*x1 + a2*x2 + a3*x3 + a4*x4
template <class OutIt, class InIt1, class InIt2, class InIt3, class InIt4, class T>
constexpr inline void scale_sum_inplace(OutIt y_first, OutIt y_last, T a1, InIt1 x1, T a2, InIt2 x2, T a3, InIt3 x3, T a4, InIt4 x4) noexcept
{
    while (y_first != y_last)
        (*y_first++) += a1 * (*x1++) + a2 * (*x2++) + a3 * (*x3++) + a4 * (*x4++);
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

} // namespace solver::it_algebra
