/// @file stepper_rk4.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "solver/detail/type_traits.hpp"
#include "solver/detail/it_algebra.hpp"

namespace solver
{

/// @brief Stepper implementing a fourth-order Runge-Kutta method.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_rk4 {
public:
    using order_type = unsigned short;
    using time_type  = Time;
    using state_type = State;
    using value_type = typename state_type::value_type;

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
        if constexpr (solver::detail::has_resize<state_type>::value) {
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
    /// @param system the system we are integrating.
    /// @param x the initial state.
    /// @param t the initial time.
    /// @param dt the step-size.
    template <class System>
    constexpr inline void do_step(System &system, state_type &x, const time_type t, const time_type dt)
    {
        const time_type dt6 = dt / static_cast<time_type>(6);
        const time_type dt3 = dt / static_cast<time_type>(3);
        const time_type dh  = dt / static_cast<time_type>(2);
        const time_type th  = t + dh;

        // m_dxdt1 = f(x, t)
        system(x, m_dxdt1, t);
        // m_x += x + dh * m_dxdt1
        detail::it_algebra::scale_sum(m_x.begin(), m_x.end(), x.begin(), dh, m_dxdt1.begin());

        // m_dxdt2 = f(m_x, t + dh)
        system(m_x, m_dxdt2, th);
        // m_x += x + dh * m_dxdt2
        detail::it_algebra::scale_sum(m_x.begin(), m_x.end(), x.begin(), dh, m_dxdt2.begin());

        // m_dxdt3 = f(m_x, t + dh)
        system(m_x, m_dxdt3, th);
        // m_x += x + dt * m_dxdt3
        detail::it_algebra::scale_sum(m_x.begin(), m_x.end(), x.begin(), dt, m_dxdt3.begin());

        // m_dxdt4 = f(m_x, t + dt)
        system(m_x, m_dxdt4, t + dt);
        // x += (dt/6)*m_dxdt1 + (dt/3)*m_dxdt2 + (dt/3)*m_dxdt3 + (dt/6)*m_dxdt4
        detail::it_algebra::scale_sum_inplace(
            x.begin(), x.end(),
            dt6, m_dxdt1.begin(),
            dt3, m_dxdt2.begin(),
            dt3, m_dxdt3.begin(),
            dt6, m_dxdt4.begin());

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

} // namespace solver