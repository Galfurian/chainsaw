/// @file stepper_improved_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/type_traits.hpp"
#include "chainsaw/detail/it_algebra.hpp"

namespace chainsaw
{

/// @brief Stepper implementing Heun's method for numerical integration (also
/// known as the Improved Euler Method).
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_improved_euler {
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
    stepper_improved_euler()
        : m_dxdt1(),
          m_dxdt2(),
          m_x(),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_improved_euler(const stepper_improved_euler &other) = delete;

    /// @brief Nope.
    stepper_improved_euler &operator=(const stepper_improved_euler &other) = delete;

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
            m_dxdt1.resize(reference.size());
            m_dxdt2.resize(reference.size());
            m_x.resize(reference.size());
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
        system(x, m_dxdt1, t);

        // Calculate the state at the next time point using Euler's method:
        //      m_x(t + dt) = x(t) + dxdt * dt;
        //
        detail::it_algebra::scale_two_sum(m_x.begin(), m_x.end(), 1., x.begin(), dt, m_dxdt1.begin());

        // Calculate the derivative at the midpoint:
        //      dxdt = system(m_x, t);
        //
        system(m_x, m_dxdt2, t + dt);

        // Increment each element of the state vector x, by half of the
        // derivative at the midpoint (dxdt) multiplied by the time step dt / 2,
        // effectively updating the state to the next time step using the
        // midpoint method:
        //      x(t + dt) = x(t) + (dt / 2) * (dxdt1 + dxdt2);
        //
        detail::it_algebra::scale_two_sum_accumulate(
            x.begin(), x.end(),
            dt * .5, m_dxdt1.begin(),
            dt * .5, m_dxdt2.begin());

        // Increase the number of steps.
        ++m_steps;
    }

private:
    /// Keeps track of state evolution.
    state_type m_dxdt1, m_dxdt2, m_x;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace chainsaw