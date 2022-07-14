/// @file stepper_adaptive_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include "stepper_euler.hpp"

#include <cmath>

namespace solver
{

template <class State, class Time>
class stepper_adaptive_euler {
public:
    using order_type_t = unsigned short;
    using time_type_t  = Time;
    using state_type_t = State;
    using value_type_t = typename State::value_type;

    stepper_adaptive_euler(value_type_t tollerance = 0.0001)
        : _stepper1(),
          _stepper2(),
          _state(),
          _tollerance(tollerance),
          _time_delta(1e-12),
          _time(.0)
    {
        // Nothing to do.
    }

    constexpr inline order_type_t order_step() const
    {
        return 0;
    }

    // Initilize the stepper.
    void initialize(const state_type_t &state, time_type_t time, time_type_t time_delta)
    {
        _state      = state;
        _time       = time;
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
    constexpr void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        _stepper1(system, x, t, dt);
    }

    template <class System>
    void do_step(System &system)
    {
        state_type_t _y0 = _state;
        // Compute values of (1) y_{n+1} = y_n + h * f(t_n, y_n).
        _stepper1.do_step(system, _y0, _time, _time_delta);
        // Compute values of (0)
        //     y_{n + 0.5} = y_n         + 0.5 * h * f(t_n, y_n)
        //     y_{n + 1}   = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
        _stepper2.do_step(system, _state, _time, _time_delta * 0.5);
        _stepper2.do_step(system, _state, _time + _time_delta * 0.5, _time_delta * 0.5);
        // Update the time.
        _time += _time_delta;
        // Update the time-delta.
        _time_delta = 0.9 * _time_delta * std::min(std::max(std::sqrt(_tollerance / (2 * __abs(_state - _y0))), 0.3), 2.0);
    }

private:
    constexpr inline auto __abs(const state_type_t &s)
    {
        return it_algebra::accumulate_abs<double>(s.begin(), s.end());
    }
    stepper_euler<State, Time> _stepper1;
    stepper_euler<State, Time> _stepper2;
    state_type_t _state;
    time_type_t _tollerance;
    time_type_t _time_delta;
    time_type_t _time;
};

} // namespace solver
