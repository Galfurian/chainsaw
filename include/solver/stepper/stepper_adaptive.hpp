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

/// Types of truncation error formulas.
enum class ErrorFormula {
    Absolute, ///< Use the absolute truncation error.
    Relative, ///< Use the relative truncation error.
    Mixed     ///< Use a mixed absolute and relative truncation error.
};

/// @brief It dynamically controlls the step-size of a stepper.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Iterations The number of iterations we are going to do while
/// integrating, higher values means more accurate results, but computationally
/// expensive.
/// @tparam Error The type of error formula we rely upon.
template <class Stepper, int Iterations = 2, ErrorFormula Error = ErrorFormula::Absolute>
class stepper_adaptive {
public:
    using stepper_type_t = Stepper;
    using order_type_t   = typename Stepper::order_type_t;
    using time_type_t    = typename Stepper::time_type_t;
    using state_type_t   = typename Stepper::state_type_t;
    using value_type_t   = typename Stepper::state_type_t::value_type;

    /// @brief Creates a new adaptive stepper.
    /// @param tollerance the tollerance we use to tweak the step-size.
    stepper_adaptive(value_type_t tollerance = 0.0001)
        : m_stepper_main(),
          m_stepper_tmp(),
          m_state(),
          m_tollerance(tollerance),
          m_time_delta(1e-12),
          m_time(.0),
          m_t_err(.0),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_adaptive(const stepper_adaptive &other) = delete;

    /// @brief Nope.
    stepper_adaptive &operator=(const stepper_adaptive &other) = delete;

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type_t order_step() const
    {
        return m_stepper_main.order_step();
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type_t &reference)
    {
        m_stepper_main.adjust_size(reference);
        m_stepper_tmp.adjust_size(reference);
    }

    /// @brief Initilizes the stepper.
    /// @param state the initial step.
    /// @param time the initial time.
    /// @param time_delta the initial step-size.
    void initialize(const state_type_t &state, time_type_t time, time_type_t time_delta)
    {
        // Initialize the state.
        m_state = state;
        // Initialize the time.
        m_time = time;
        // Initialize the step-size.
        m_time_delta = time_delta;
    }

    /// @brief Returns a copy of the current state.
    /// @return a copy of the state.
    inline state_type_t current_state() const
    {
        return m_state;
    }

    /// @brief Returns a copy of the current step-size.
    /// @return a copy of the step-size.
    inline time_type_t current_time_step() const
    {
        return m_time_delta;
    }

    /// @brief Returns a copy of the current time.
    /// @return a copy of the time.
    inline time_type_t current_time() const
    {
        return m_time;
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
    constexpr inline void do_step(System &system, state_type_t &x, time_type_t t, time_type_t dt) noexcept
    {
        // Call the stepper.
        m_stepper_main(system, x, t, dt);
        // Increase the number of steps.
        ++m_steps;
    }

    /// @brief Integrates on step.
    /// @param system the system we are integrating.
    template <class System>
    constexpr inline void do_step(System &system)
    {
        state_type_t m_y0(m_state);

        // Compute values of (1) y_{n+1} = y_n + h * f(t_n, y_n).
        m_stepper_main.do_step(system, m_y0, m_time, m_time_delta);

        // Compute values of (0)
        //     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
        //     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
        if constexpr (Iterations <= 2) {
            const time_type_t dh = m_time_delta * .5;
            m_stepper_tmp.do_step(system, m_state, m_time + dh, dh);
            m_stepper_tmp.do_step(system, m_state, m_time + m_time_delta, dh);
        } else {
            constexpr time_type_t hs = 1. / Iterations;
            const time_type_t dh     = m_time_delta * hs;
            for (int i = 0; i < Iterations; ++i)
                m_stepper_tmp.do_step(system, m_state, m_time + dh * i, dh);
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
    /// The main stepper.
    stepper_type_t m_stepper_main;
    /// A temporary stepper we use to tune the main stepper.
    stepper_type_t m_stepper_tmp;
    /// The current step.
    state_type_t m_state;
    /// The tollerance value we use to tune the step-size.
    time_type_t m_tollerance;
    /// The current step-size, it is not fixed.
    time_type_t m_time_delta;
    /// The current time.
    time_type_t m_time;
    /// Holds the error between the main stepper and the temporary stepper.
    value_type_t m_t_err;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace solver
