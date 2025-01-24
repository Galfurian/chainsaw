/// @file less_with_sign.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

namespace numint::detail
{

/// @brief Returns (t1 < t2) if (dt > 0) and (t1 > t2) if (dt < 0) with epsilon accuracy.
/// @param t1
/// @param t2
/// @param dt
/// @return (t1 < t2) if (dt > 0) and (t1 > t2) if (dt < 0) with epsilon accuracy.
template <typename T>
inline bool less_with_sign(T t1, T t2, T dt)
{
    return (dt > 0) ? t1 < t2 : t2 < t1;
}

/// @brief Returns (t1 <= t2) if (dt > 0) and (t1 >= t2) if (dt < 0) with epsilon accuracy.
/// @param t1
/// @param t2
/// @param dt
/// @return (t1 <= t2) if (dt > 0) and (t1 >= t2) if (dt < 0) with epsilon accuracy.
template <typename T>
inline bool less_eq_with_sign(T t1, T t2, T dt)
{
    return (dt > 0) ? t1 <= t2 : t2 <= t1;
}

} // namespace numint::detail
