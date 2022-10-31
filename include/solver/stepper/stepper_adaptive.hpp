/// @file stepper_adaptive.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "solver/detail/type_traits.hpp"
#include "solver/detail/it_algebra.hpp"

#include <cmath>

namespace solver
{

enum class ErrorFormula {
    Absolute,
    Relative,
    Mixed
};

template <class State, class Time, class Stepper, int Iterations = 2, ErrorFormula Error = ErrorFormula::Absolute>
class stepper_adaptive {
public:
    using order_type_t   = unsigned short;
    using time_type_t    = Time;
    using state_type_t   = State;
    using stepper_type_t = Stepper;
    using value_type_t   = typename State::value_type;

    stepper_adaptive(value_type_t tollerance = 0.0001)
        : m_stepper1(),
          m_stepper2(),
          m_state(),
          m_tollerance(tollerance),
          m_time_delta(1e-12),
          m_time(.0),
          m_t_err(.0),
          m_steps()
    {
        // Nothing to do.
    }

    stepper_adaptive(const stepper_adaptive &other) = delete;

    stepper_adaptive &operator=(const stepper_adaptive &other) = delete;

    constexpr inline order_type_t order_step() const
    {
        return 0;
    }

    void adjust_size(state_type_t &x)
    {
        m_stepper1.adjust_size(x);
        m_stepper2.adjust_size(x);
    }

    // Initilize the stepper.
    void initialize(const state_type_t &state, time_type_t time, time_type_t time_delta)
    {
        // Initialize the state.
        m_state = state;
        // Initialize the time.
        m_time = time;
        // Initialize the step size.
        m_time_delta = time_delta;
    }

    inline state_type_t current_state() const
    {
        return m_state;
    }

    inline time_type_t current_time_step() const
    {
        return m_time_delta;
    }

    inline time_type_t current_time() const
    {
        return m_time;
    }

    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr inline void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        m_stepper1(system, x, t, dt);

        // Increase the number of steps.
        ++m_steps;
    }

    template <class System>
    constexpr inline void do_step(System &system)
    {
        state_type_t m_y0 = m_state;

        // Compute values of (1) y_{n+1} = y_n + h * f(t_n, y_n).
        m_stepper1.do_step(system, m_y0, m_time, m_time_delta);

        // Compute values of (0)
        //     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
        //     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
        if constexpr (Iterations <= 2) {
            const Time dh = m_time_delta * .5;
            m_stepper2.do_step(system, m_state, m_time + dh, dh);
            m_stepper2.do_step(system, m_state, m_time + m_time_delta, dh);
        } else {
            constexpr Time hs = 1. / Iterations;
            const Time dh     = m_time_delta * hs;
            for (int i = 0; i < Iterations; ++i)
                m_stepper2.do_step(system, m_state, m_time + dh * i, dh);
        }

        // Update the time.
        m_time += m_time_delta;

        // Calculate truncation error.
        if constexpr (Error == ErrorFormula::Absolute) {
            // Use absolute truncation error.
            m_t_err = detail::it_algebra::max_abs_diff<value_type_t>(m_state.begin(), m_state.end(), m_y0.begin(), m_y0.end());
        } else if constexpr (Error == ErrorFormula::Relative) {
            // Use relative truncation error.
            m_t_err = detail::it_algebra::max_rel_diff<value_type_t>(m_state.begin(), m_state.end(), m_y0.begin(), m_y0.end());
        } else {
            // Use mixed truncation error.
            m_t_err = detail::it_algebra::max_comb_diff<value_type_t>(m_state.begin(), m_state.end(), m_y0.begin(), m_y0.end());
        }

        // Update the time-delta.
        m_time_delta *= 0.9 * std::min(std::max(std::pow(m_tollerance / (2 * m_t_err), 0.2), 0.3), 2.);

        // Increase the number of steps.
        ++m_steps;
    }

private:
    stepper_type_t m_stepper1;
    stepper_type_t m_stepper2;
    state_type_t m_state;
    time_type_t m_tollerance;
    time_type_t m_time_delta;
    time_type_t m_time;
    value_type_t m_t_err;
    unsigned long m_steps;
};

} // namespace solver
