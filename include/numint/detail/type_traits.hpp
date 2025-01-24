/// @file type_traits.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Helper functions for dealing with template types.

#pragma once

#include <type_traits>

namespace numint::detail
{

/// @brief Checks if a type has a resize method.
/// @tparam T The type to check.
template <typename T, typename = void>
struct has_resize : std::false_type {
};

/// @brief Checks if a type has a resize method. 
/// @tparam T The type to check.
template <typename T>
struct has_resize<T, std::void_t<decltype(std::declval<T>().resize(1))>> : std::true_type {
};

/// @brief Helper variable template to check if a type has a resize method.
/// @tparam T The type to check.
template <typename T>
constexpr inline bool has_resize_v = has_resize<T>::value;

}