/// @file type_traits.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Helper functions for dealing with template types.

#pragma once

#include <type_traits>

template <typename T, typename = int>
struct has_resize : std::false_type {
};

template <typename T>
struct has_resize<T, decltype((void)std::declval<T>().resize(1), 0)> : std::true_type {
};

template <typename T>
constexpr inline bool has_resize_v = has_resize<T>::value;
