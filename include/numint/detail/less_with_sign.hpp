/// @file less_with_sign.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

namespace numint::detail
{

/// @brief Compares two values considering the sign of a direction parameter.
///
/// @details Returns (value1 < value2) if (direction > 0) and (value1 > value2)
/// if (direction < 0), ensuring a consistent comparison with the given direction.
///
/// @param value1 The first value to compare.
/// @param value2 The second value to compare.
/// @param direction The direction parameter determining the comparison order.
/// @return True if the condition is met, false otherwise.
template <typename T>
inline auto less_with_sign(T value1, T value2, T direction) -> bool
{
    return (direction > 0) ? value1 < value2 : value2 < value1;
}

/// @brief Compares two values with equality consideration, based on the direction parameter.
///
/// @details Returns (value1 <= value2) if (direction > 0) and (value1 >= value2)
/// if (direction < 0), ensuring a consistent comparison with the given direction.
///
/// @param value1 The first value to compare.
/// @param value2 The second value to compare.
/// @param direction The direction parameter determining the comparison order.
/// @return True if the condition is met, false otherwise.
template <typename T>
inline auto less_eq_with_sign(T value1, T value2, T direction) -> bool
{
    return (direction > 0) ? value1 <= value2 : value2 <= value1;
}

} // namespace numint::detail
