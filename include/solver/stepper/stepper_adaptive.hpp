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
    using stepper_type = Stepper;
    using order_type   = typename Stepper::order_type;
    using time_type    = typename Stepper::time_type;
    using state_type   = typename Stepper::state_type;
    using value_type   = typename Stepper::state_type::value_type;

    /// @brief Creates a new adaptive stepper.
    /// @param tollerance the tollerance we use to tweak the step-size.
    stepper_adaptive()
        : _stepper_main(),
          _stepper_tuner(),
          _tollerance(0.0001),
          _min_delta(1e-12),
          _max_delta(1),
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

    constexpr inline void set_tollerance(value_type tollerance)
    {
        _tollerance = tollerance;
    }

    constexpr inline void set_min_delta(value_type min_delta)
    {
        _min_delta = min_delta;
    }

    constexpr inline void set_max_delta(value_type max_delta)
    {
        _max_delta = max_delta;
    }

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type order_step() const
    {
        return _stepper_main.order_step();
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type &reference)
    {
        _stepper_main.adjust_size(reference);
        _stepper_tuner.adjust_size(reference);
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr inline auto steps() const
    {
        return _steps;
    }

    /// @brief Integrates on step.
    /// @param system the system we are integrating.
    /// @details
    /// Compute values of (0)
    ///     y_{n + 1}   = y_n + h * f(t_n, y_n)
    /// Compute values of (1)
    ///     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
    ///     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
    ///
    template <class System>
    constexpr inline void do_step(System &system, state_type &x, time_type &t, time_type &dt)
    {
        state_type y(x);
        // Compute values of (0).
        _stepper_main.do_step(system, y, t, dt);
        // Compute values of (1).
        if constexpr (Iterations <= 2) {
            const time_type dh = dt * .5;
            _stepper_tuner.do_step(system, x, t + dh, dh);
            _stepper_tuner.do_step(system, x, t + dt, dh);
        } else {
            constexpr time_type hs = 1. / Iterations;
            const time_type dh     = dt * hs;
            for (int i = 0; i < Iterations; ++i)
                _stepper_tuner.do_step(system, x, t + dh * i, dh);
        }
        // Update the time.
        t += dt;
        // Calculate truncation error.
        if constexpr (Error == ErrorFormula::Absolute) {
            // Get absolute truncation error.
            _t_err = _t_err_abs = detail::it_algebra::max_abs_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            dt *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err_abs), 0.2), 0.3), 2.);
        } else if constexpr (Error == ErrorFormula::Relative) {
            // Get relative truncation error.
            _t_err = _t_err_rel = detail::it_algebra::max_rel_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            dt *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err_rel), 0.2), 0.3), 2.);
        } else {
#if 1
            // Get mixed truncation error.
            _t_err = detail::it_algebra::max_comb_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            dt *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err), 0.2), 0.3), 2.);
#else
            // Get absolute truncation error.
            _t_err_abs = detail::it_algebra::max_abs_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Get relative truncation error.
            _t_err_rel = detail::it_algebra::max_rel_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Get mixed truncation error.
            _t_err = _t_err_abs + _t_err_rel;
            // Update the time-delta.
            dt *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err), 0.2), 0.2), 4.);
#endif
        }
        // Check boundaries.
        dt = std::min(std::max(dt, _min_delta), _max_delta);
        // Increase the number of steps.
        ++_steps;
    }

private:
    /// The main stepper.
    stepper_type _stepper_main;
    /// A temporary stepper we use to tune the main stepper.
    stepper_type _stepper_tuner;
    /// The tollerance value we use to tune the step-size.
    time_type _tollerance;
    /// The minimum step-size.
    time_type _min_delta;
    /// The maximum step-size.
    time_type _max_delta;
    /// Holds the error between the main stepper and the temporary stepper.
    value_type _t_err;
    /// Holds the absolute error between the main stepper and the temporary stepper.
    value_type _t_err_abs;
    /// Holds the relative error between the main stepper and the temporary stepper.
    value_type _t_err_rel;
    /// The number of steps of integration.
    unsigned long _steps;
};

} // namespace solver
