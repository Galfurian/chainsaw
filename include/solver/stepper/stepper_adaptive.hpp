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
        : _stepper1(),
          _stepper2(),
          _state(),
          _tollerance(tollerance),
          _time_delta(1e-12),
          _time(.0),
          _t_err(.0)
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
        _stepper1.adjust_size(x);
        _stepper2.adjust_size(x);
    }

    // Initilize the stepper.
    void initialize(const state_type_t &state, time_type_t time, time_type_t time_delta)
    {
        // Initialize the state.
        _state = state;
        // Initialize the time.
        _time = time;
        // Initialize the step size.
        _time_delta = time_delta;
    }

    inline state_type_t current_state() const
    {
        return _state;
    }

    inline time_type_t current_time_step() const
    {
        return _time_delta;
    }

    inline time_type_t current_time() const
    {
        return _time;
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr inline void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        _stepper1(system, x, t, dt);
    }

    template <class System>
    constexpr inline void do_step(System &system)
    {
        state_type_t _y0 = _state;

        // Compute values of (1) y_{n+1} = y_n + h * f(t_n, y_n).
        _stepper1.do_step(system, _y0, _time, _time_delta);

        // Compute values of (0)
        //     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
        //     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
        if constexpr (Iterations <= 2) {
            const Time dh = _time_delta * .5;
            _stepper2.do_step(system, _state, _time + dh, dh);
            _stepper2.do_step(system, _state, _time + _time_delta, dh);
        } else {
            constexpr Time hs = 1. / Iterations;
            const Time dh     = _time_delta * hs;
            for (int i = 0; i < Iterations; ++i)
                _stepper2.do_step(system, _state, _time + dh * i, dh);
        }

        // Update the time.
        _time += _time_delta;

        // Calculate truncation error.
        if constexpr (Error == ErrorFormula::Absolute) {
            // Use absolute truncation error.
            _t_err = detail::it_algebra::max_abs_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
        } else if constexpr (Error == ErrorFormula::Relative) {
            // Use relative truncation error.
            _t_err = detail::it_algebra::max_rel_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
        } else {
            // Use mixed truncation error.
            _t_err = detail::it_algebra::max_comb_diff<value_type_t>(_state.begin(), _state.end(), _y0.begin(), _y0.end());
        }

        // Update the time-delta.
        _time_delta *= 0.9 * std::min(std::max(std::pow(_tollerance / (2 * _t_err), 0.2), 0.3), 2.);
    }

private:
    stepper_type_t _stepper1;
    stepper_type_t _stepper2;
    state_type_t _state;
    time_type_t _tollerance;
    time_type_t _time_delta;
    time_type_t _time;
    value_type_t _t_err;
    unsigned _iterations;
};

} // namespace solver
