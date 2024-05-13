/// @file stepper_trapezoidal.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/type_traits.hpp"
#include "chainsaw/detail/it_algebra.hpp"

namespace chainsaw
{

/// @brief Stepper implementing the trapezoidal method.
/// @details Approximates the area under the curve by dividing the interval into trapezoids.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_trapezoidal {
public:
    /// @brief Type used for the order of the stepper.
    using order_type = unsigned short;
    /// @brief Type used to keep track of time.
    using time_type = Time;
    /// @brief The state vector.
    using state_type = State;
    /// @brief Type of value contained in the state vector.
    using value_type = typename state_type::value_type;
    /// @brief Determines if this is an adaptive stepper or not.
    static constexpr bool is_adaptive_stepper = false;

    /// @brief Creates a new stepper.
    stepper_trapezoidal()
        : m_dxdt_start(),
          m_dxdt_end(),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_trapezoidal(const stepper_trapezoidal &other) = delete;

    /// @brief Nope.
    stepper_trapezoidal &operator=(const stepper_trapezoidal &other) = delete;

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type order_step() const
    {
        return 1;
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type &reference)
    {
        if constexpr (chainsaw::detail::has_resize<state_type>::value) {
            m_dxdt_start.resize(reference.size());
            m_dxdt_end.resize(reference.size());
        }
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Perform a single integration step using Euler's method.
    /// @tparam System The type of the system representing the differential equations.
    /// @param system the system we are integrating.
    /// @param x the initial state.
    /// @param t the initial time.
    /// @param dt the step-size.
    template <class System>
    constexpr void do_step(System &system, state_type &x, const time_type t, const time_type dt) noexcept
    {
        // Calculate the derivative at the start point.
        //
        system(x, m_dxdt_start, t);

        // Calculate the derivative at the end point.
        //
        system(x, m_dxdt_end, t + dt);

        // Update the state vector using Euler's method:
        //      x(t + dt) = x(t) + (0.5 * dt * dxdt_start) + (0.5 * dt * dxdt_end)
        detail::it_algebra::scale_two_sum_accumulate(
            x.begin(), x.end(),
            0.5 * dt, m_dxdt_start.begin(),
            0.5 * dt, m_dxdt_end.begin());

        // Increment the number of integration steps.
        ++m_steps;
    }

private:
    /// Keeps track of state evolution.
    state_type m_dxdt_start, m_dxdt_end;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace chainsaw