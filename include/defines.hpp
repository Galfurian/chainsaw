/// @file defines.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-04-13

#pragma once

#include <iostream>
#include <array>

/// @brief Time itself.
using Time = double;

/// @brief Returns the sign of the input value (+1 or -1).
template <typename T>
constexpr inline int sign(const T &value)
{
    return (value > 0) - (value < 0);
}

template <typename T, size_t size>
std::ostream &operator<<(std::ostream &lhs, const std::array<T, size> &rhs)
{
    lhs << "[ ";
    for (size_t i = 0; i < size; ++i)
        std::cout << rhs[i] << " ";
    lhs << "]";
    return lhs;
}