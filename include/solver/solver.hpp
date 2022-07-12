/// @file solver.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include "it_algebra.hpp"
#include "stepper_euler.hpp"
#include "stepper_rk4.hpp"
#include "stepper_adaptive.hpp"
#include "less_with_sign.hpp"

namespace solver
{

namespace detail
{

template <class Stepper, class System, class Observer>
inline constexpr void integrate_one_step_const(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t time,
    typename Stepper::time_type_t time_delta) noexcept
{
    stepper.do_step(system, state, time, time_delta);
    observer(state, time);
}

} // namespace detail

template <class Stepper, class System, class Observer>
inline constexpr unsigned integrate_const(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta) noexcept
{
    unsigned iteration = 0;
    observer(state, start_time);
    while (start_time < end_time) {
        detail::integrate_one_step_const(stepper, observer, system, state, start_time, time_delta);
        start_time += time_delta;
        ++iteration;
    }
    return iteration;
}

template <class Stepper, class System, class Observer>
inline constexpr unsigned integrate_adaptive_simple(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta) noexcept
{
    unsigned iteration = integrate_const(stepper, observer, system, state, start_time, end_time, time_delta);
    // Make a last step to end exactly at end_time.
    if (less_with_sign(start_time + time_delta * iteration, end_time, time_delta)) {
        detail::integrate_one_step_const(stepper, observer, system, state, start_time, time_delta);
        ++iteration;
    }
    return iteration;
}

template <class Stepper, class System, class Observer>
unsigned integrate_adaptive(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta)
{
    unsigned iteration = 0;
    // Initilize the stepper.
    stepper.initialize(state, start_time, time_delta);
    while (less_with_sign(stepper.current_time(), end_time, stepper.current_time_step())) {
        // Make sure we don't go beyond the end_time.
        while (less_eq_with_sign(static_cast<Stepper::time_type_t>(stepper.current_time() + stepper.current_time_step()), end_time, stepper.current_time_step())) {
            stepper.do_step(system);
            observer(stepper.current_state(), stepper.current_time());
            ++iteration;
        }
        // Calculate time step to arrive exactly at end time.
        stepper.initialize(stepper.current_state(), stepper.current_time(), static_cast<Stepper::time_type_t>(end_time - stepper.current_time()));
    }
    observer(stepper.current_state(), stepper.current_time());
    // Overwrite state with the final point.
    state = stepper.current_state();
    return iteration;
}

} // namespace solver