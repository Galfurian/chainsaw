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
    stepper_adaptive()
        : _stepper_main(),
          _stepper_tmp(),
          _state(),
          _tollerance(0.0001),
          _time_delta(1e-12),
          _min_delta(1e-12),
          _max_delta(1),
          _time(.0),
          _t_err(.0),
          _t_err_abs(.0),
          _t_err_rel(.0),
          _steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_adaptive(const stepper_adaptive &other) = delete;

    /// @brief Nope.
    stepper_adaptive &operator=(const stepper_adaptive &other) = delete;

    constexpr inline void set_tollerance(value_type_t tollerance)
    {
        _tollerance = tollerance;
    }

    constexpr inline void set_min_delta(value_type_t min_delta)
    {
        _min_delta = min_delta;
    }

    constexpr inline void set_max_delta(value_type_t max_delta)
    {
        _max_delta = max_delta;
    }

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type_t order_step() const
    {
        return _stepper_main.order_step();
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type_t &reference)
    {
        _stepper_main.adjust_size(reference);
        _stepper_tmp.adjust_size(reference);
    }

    /// @brief Initilizes the stepper.
    /// @param state the initial step.
    /// @param time the initial time.
    /// @param time_delta the initial step-size.
    void initialize(const state_type_t &state, time_type_t time, time_type_t time_delta)
    {
        // Initialize the state.
        _state = state;
        // Initialize the time.
        _time = time;
        // Initialize the step-size.
        _time_delta = time_delta;
    }

    /// @brief Returns a copy of the current state.
    /// @return a copy of the state.
    inline state_type_t current_state() const
    {
        return _state;
    }

    /// @brief Returns a copy of the current step-size.
    /// @return a copy of the step-size.
    inline time_type_t current_time_step() const
    {
        return _time_delta;
    }

    /// @brief Returns a copy of the current time.
    /// @return a copy of the time.
    inline time_type_t current_time() const
    {
        return _time;
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr inline auto steps() const
    {
        return _steps;
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
        _stepper_main(system, x, t, dt);
        // Increase the number of steps.
        ++_steps;
    }

    /// @brief Integrates on step.
    /// @param system the system we are integrating.
    template <class System>
    constexpr inline void do_step(System &system)
    {
        state_type_t _y0(_state);

        // Compute values of (1) y_{n+1} = y_n + h * f(t_n, y_n).
        _stepper_main.do_step(system, _y0, _time, _time_delta);

        // Compute values of (0)
        //     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
        //     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
        if constexpr (Iterations <= 2) {
            const time_type_t dh = _time_delta * .5;
            _stepper_tmp.do_step(system, _state, _time + dh, dh);
            _stepper_tmp.do_step(system, _state, _time + _time_delta, dh);
        } else {
            constexpr time_type_t hs = 1. / Iterations;
            const time_type_t dh     = _time_delta * hs;
            for (int i = 0; i < Iterations; ++i)
                _stepper_tmp.do_step(system, _state, _time + dh * i, dh);
        }

        // Update the time.
        _time += _time_delta;

        // Calculate truncation error.
        if constexpr (Error == ErrorFormula::Absolute) {
            // Get absolute truncation error.
            _t_err = _t_err_abs = detail::it_algebra::max_abs_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
            // Update the time-delta.
            _time_delta *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err_abs), 0.2), 0.3), 2.);
        } else if constexpr (Error == ErrorFormula::Relative) {
            // Get relative truncation error.
            _t_err = _t_err_rel = detail::it_algebra::max_rel_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
            // Update the time-delta.
            _time_delta *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err_rel), 0.2), 0.3), 2.);
        } else {
#if 1
            // Get mixed truncation error.
            _t_err = detail::it_algebra::max_comb_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
            // Update the time-delta.
            _time_delta *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err), 0.2), 0.3), 2.);
#else
            // Get absolute truncation error.
            _t_err_abs = detail::it_algebra::max_abs_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
            // Get relative truncation error.
            _t_err_rel = detail::it_algebra::max_rel_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
            // Get mixed truncation error.
            _t_err = _t_err_abs + _t_err_rel;
            // Update the time-delta.
            _time_delta *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err), 0.2), 0.2), 4.);
#endif
        }
        // Check boundaries.
        _time_delta = std::min(std::max(_time_delta, _min_delta), _max_delta);

        // Increase the number of steps.
        ++_steps;
    }

private:
    /// The main stepper.
    stepper_type_t _stepper_main;
    /// A temporary stepper we use to tune the main stepper.
    stepper_type_t _stepper_tmp;
    /// The current step.
    state_type_t _state;
    /// The tollerance value we use to tune the step-size.
    time_type_t _tollerance;
    /// The current step-size, it is not fixed.
    time_type_t _time_delta;
    /// The minimum step-size.
    time_type_t _min_delta;
    /// The maximum step-size.
    time_type_t _max_delta;
    /// The current time.
    time_type_t _time;
    /// Holds the error between the main stepper and the temporary stepper.
    value_type_t _t_err;
    /// Holds the absolute error between the main stepper and the temporary stepper.
    value_type_t _t_err_abs;
    /// Holds the relative error between the main stepper and the temporary stepper.
    value_type_t _t_err_rel;
    /// The number of steps of integration.
    unsigned long _steps;
};

} // namespace solver
