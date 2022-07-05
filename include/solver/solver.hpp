/// @file solver.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include "it_algebra.hpp"
#include "stepper_euler.hpp"
#include "stepper_rk4.hpp"

namespace solver
{

template <class Stepper, class System, class Observer>
constexpr void integrate_one_step(
    Stepper &stepper,
    System &system,
    typename Stepper::container_type_t &state,
    typename Stepper::time_type_t &time,
    typename Stepper::time_type_t time_delta,
    Observer &observer) noexcept
{
    stepper.do_step(system, state, time, time_delta);
    observer(state, time);
    time += time_delta;
}

template <class Stepper, class System, class Observer>
constexpr unsigned integrate_const(
    Stepper &stepper,
    System &system,
    typename Stepper::container_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta,
    Observer &observer) noexcept
{
    unsigned iteration = 0;
    observer(state, start_time);
    while (start_time < end_time) {
        integrate_one_step(stepper, system, state, start_time, time_delta, observer);
        ++iteration;
    }
    return iteration;
}

} // namespace solver