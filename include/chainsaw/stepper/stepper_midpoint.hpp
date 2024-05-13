/// @file stepper_midpoint.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/type_traits.hpp"
#include "chainsaw/detail/it_algebra.hpp"

namespace chainsaw
{

/// @brief Stepper implementing Midpoint Method (Modified Euler's Method).
/// @tparam Stepper The stepper we rely upon.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_midpoint {
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
    stepper_midpoint()
        : m_dxdt(),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_midpoint(const stepper_midpoint &other) = delete;

    /// @brief Nope.
    stepper_midpoint &operator=(const stepper_midpoint &other) = delete;

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
            m_dxdt.resize(reference.size());
        }
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Integrates on step.
    /// @param system the system we are integrating.
    /// @param x the initial state.
    /// @param t the initial time.
    /// @param dt the step-size.
    template <class System>
    constexpr void do_step(System &system, state_type &x, const time_type t, const time_type dt) noexcept
    {
        // Calculate the derivative at the initial point:
        //      dxdt = system(x, t);
        // 
        system(x, m_dxdt, t);

        // Increment each element of the state vector x, by half of the
        // derivative at the initial point (dxdt) multiplied by the time step dt / 2,
        // effectively updating the state to the midpoint:
        //      x(t + (dt / 2)) = x(t) + dxdt * (dt / 2);
        // 
        detail::it_algebra::scale_accumulate(x.begin(), x.end(), m_dxdt.begin(), dt / 2.);
        
        // Calculate the derivative at the midpoint:
        //      dxdt = system(x, t);
        // 
        system(x, m_dxdt, t + (dt / 2.));

        // Increment each element of the state vector x, by half of the
        // derivative at the midpoint (dxdt) multiplied by the time step dt / 2,
        // effectively updating the state to the next time step using the
        // midpoint method:
        //      x(t + dt) = x(t) + dxdt * (dt / 2);
        // 
        detail::it_algebra::scale_accumulate(x.begin(), x.end(), m_dxdt.begin(), dt / 2.);
        
        // Increase the number of steps.
        ++m_steps;
    }

private:
    /// Keeps track of state evolution.
    state_type m_dxdt;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace chainsaw