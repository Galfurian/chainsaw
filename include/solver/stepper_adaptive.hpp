/// @file stepper_adaptive.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include "defines.hpp"

#include <exception>

namespace solver
{

template <class State, class Time>
class stepper_adaptive {
public:
    using order_type_t = unsigned short;
    using time_type_t  = Time;
    using state_type_t = State;
    using value_type_t = typename State::value_type;

    stepper_adaptive(time_type_t time_delta_min, time_type_t time_delta_max, int integration_steps = 100, value_type_t max_rate_of_change = 0.001)
        : _state(),
          _dxdt(),
          _time_delta(time_delta_min),
          _time_delta_min(time_delta_min),
          _time_delta_max(time_delta_max),
          _time(.0),
          _max_rate_of_change(max_rate_of_change),
          _goodsteps(0),
          _integration_steps(integration_steps),
          _stepcount(0)
    {
        if (time_delta_min > time_delta_max)
            throw std::runtime_error("Minimum time_delta should be lower than maximum");
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

    inline int get_steps() const
    {
        return _stepcount;
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        system(x, _dxdt, t);
        it_algebra::increment(x.begin(), x.end(), _dxdt.begin(), dt);
    }

    template <class System>
    void do_step(System &system)
    {
        state_type_t dxdt1, dxdt2;

        ++_stepcount;

        // Simulate the system.
        system(_state, dxdt1, _time);
        system(_state + dxdt1 * (_time_delta * .5), dxdt2, _time + (_time_delta * .5));

        auto abs_e1 = abs(_state + dxdt1 * _time_delta);
        auto abs_e2 = abs(_state + dxdt1 * (_time_delta * .5) + dxdt2 * (_time_delta * .5));

        auto rel_change = std::abs(abs_e1 - abs_e2) / (std::abs(abs_e1 + abs_e2) * .5);

        if (rel_change >= _max_rate_of_change) {
            _time_delta = std::max(_time_delta * .5, _time_delta_min);
            _goodsteps = 0;
        } else {
            ++_goodsteps;

            if (_goodsteps >= _integration_steps){
                _time_delta = std::min(_time_delta * 1.1, _time_delta_max);
                _goodsteps = 0;
            }
        }

        _state += dxdt1 * _time_delta;
        _time += _time_delta;

    }

private:
    state_type_t _state;
    state_type_t _dxdt;
    time_type_t _time_delta;
    time_type_t _time_delta_min;
    time_type_t _time_delta_max;
    time_type_t _time;
    value_type_t _max_rate_of_change;
    int _goodsteps;
    int _integration_steps;
    int _stepcount;
};

} // namespace solver
