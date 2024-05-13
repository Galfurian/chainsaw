/// @file stepper_rk4.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/type_traits.hpp"
#include "chainsaw/detail/it_algebra.hpp"

namespace chainsaw
{

/// @brief Stepper implementing a fourth-order Runge-Kutta method.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_rk4 {
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
    stepper_rk4()
        : m_dxdt1(),
          m_dxdt2(),
          m_dxdt3(),
          m_dxdt4(),
          m_x(),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_rk4(const stepper_rk4 &other) = delete;

    /// @brief Nope.
    stepper_rk4 &operator=(const stepper_rk4 &other) = delete;

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type order_step() const
    {
        return 4;
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    constexpr inline void adjust_size(const state_type &reference)
    {
        if constexpr (chainsaw::detail::has_resize<state_type>::value) {
            m_dxdt1.resize(reference.size());
            m_dxdt2.resize(reference.size());
            m_dxdt3.resize(reference.size());
            m_dxdt4.resize(reference.size());
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
    /// @tparam System The type of the system representing the differential equations.
    /// @param system the system we are integrating.
    /// @param x the initial state.
    /// @param t the initial time.
    /// @param dt the step-size.
    template <class System>
    constexpr inline void do_step(System &system, state_type &x, const time_type t, const time_type dt)
    {
        // Here is the idea:
        //  - m_dxdt1 : Slope at the beginning of the interval
        //  - m_dxdt2 : Slope at the midpoint of the interval
        //  - m_dxdt3 : Another slope at the midpoint of the interval
        //  - m_dxdt4 : Slope at the end of the interval

        // Step 1: Calculate the slope at the beginning of the interval (m_dxdt1):
        //      m_dxdt1 = f(x, t);
        //
        system(x, m_dxdt1, t);

        // Update temporary state using the slope at the beginning and move halfway forward:
        //      m_x(t + dt * 0.5) = x(t) + m_dxdt1 * dt * 0.5;
        //
        detail::it_algebra::scale_two_sum(m_x.begin(), m_x.end(), 1., x.begin(), 0.5 * dt, m_dxdt1.begin());

        // Step 2: Calculate the slope at the midpoint of the interval (m_dxdt2):
        //      m_dxdt2 = f(m_x, t + 0.5 * dt);
        //
        system(m_x, m_dxdt2, t + 0.5 * dt);

        // Update temporary state using the slope at the midpoint and move halfway forward again:
        //      m_x(t + dt * 0.5) = x(t) + m_dxdt2 * dt * 0.5;
        //
        detail::it_algebra::scale_two_sum(m_x.begin(), m_x.end(), 1., x.begin(), 0.5 * dt, m_dxdt2.begin());

        // Step 3: Calculate another slope at the midpoint of the interval (m_dxdt3):
        //      m_dxdt3 = f(m_x, t + 0.5 * dt);
        //
        system(m_x, m_dxdt3, t + 0.5 * dt);

        // Update temporary state using the slope at the midpoint and move to the end of the interval:
        //      m_x(t + dt) = x(t) + m_dxdt3 * dt;
        //
        detail::it_algebra::scale_two_sum(m_x.begin(), m_x.end(), 1., x.begin(), dt, m_dxdt3.begin());

        // Step 4: Calculate the slope at the end of the interval (m_dxdt4):
        //      m_dxdt4 = f(m_x, t + dt);
        system(m_x, m_dxdt4, t + dt);

        // Update each component of the state vector using the weighted average
        // of the slopes: m_dxdt1, m_dxdt2, m_dxdt3, and m_dxdt4.
        detail::it_algebra::scale_four_sum_accumulate(
            x.begin(), x.end(),
            dt * (1. / 6.), m_dxdt1.begin(),
            dt * (2. / 6.), m_dxdt2.begin(),
            dt * (2. / 6.), m_dxdt3.begin(),
            dt * (1. / 6.), m_dxdt4.begin());

        // Increase the number of steps.
        ++m_steps;
    }

private:
    /// First support vector.
    state_type m_dxdt1;
    /// Second support vector.
    state_type m_dxdt2;
    /// Third support vector.
    state_type m_dxdt3;
    /// Fourth support vector.
    state_type m_dxdt4;
    /// Temporary state vector.
    state_type m_x;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace chainsaw