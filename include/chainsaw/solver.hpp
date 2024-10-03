/// @file solver.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/less_with_sign.hpp"
#include "chainsaw/detail/it_algebra.hpp"
#include "chainsaw/detail/type_traits.hpp"

namespace chainsaw
{

namespace detail
{

/// @brief Calls the integration stepper with the given system.
/// @param stepper the stepper we are using.
/// @param observer the observer we need to call after executing one integration step.
/// @param system the system we are integrating.
/// @param state the state of the system we need to evolve for one step.
/// @param time the initial time.
/// @param time_delta the integration step size.
template <class Stepper, class System, class Observer>
constexpr inline void integrate_one_step(
    Stepper &stepper,
    Observer &observer,
    System &&system,
    typename Stepper::state_type &state,
    const typename Stepper::time_type time,
    const typename Stepper::time_type time_delta) noexcept
{
    // Perform one integration step.
    stepper.do_step(system, state, time, time_delta);
    // Call the observer.
    observer(state, time);
}

} // namespace detail

/// @brief Integrates the system between the start and end time, with the given stepper.
/// @param stepper the stepper we are using.
/// @param observer the observer we need to call after executing one integration step.
/// @param system the system we are integrating.
/// @param state the initial state of the system.
/// @param start_time the start time.
/// @param end_time the final time.
/// @param time_delta the fixed integration step.
/// @return the number of steps it took to perform the integration.
template <class Stepper, class System, class Observer>
constexpr inline auto integrate_fixed(
    Stepper &stepper,
    Observer &observer,
    System &&system,
    typename Stepper::state_type &state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type time_delta) noexcept
{
    using state_type = typename Stepper::state_type;
    // Check if the state vector can (and should) be resized.
    if constexpr (chainsaw::detail::has_resize_v<state_type>) {
        stepper.adjust_size(state);
    }
    // Call the observer at the beginning.
    observer(state, start_time);
    // Run until the time reaches the `end_time`.
    while (start_time < end_time) {
        // Integrate one step.
        detail::integrate_one_step(stepper, observer, system, state, start_time, time_delta);
        // Advance time.
        start_time += time_delta;
    }
    // Return the number of steps it took to integrate.
    return stepper.steps();
}

/// @brief Integrates the system between the start and end time, with the given stepper.
/// @param stepper the stepper we are using.
/// @param observer the observer we need to call after executing one integration step.
/// @param system the system we are integrating.
/// @param state the initial state of the system.
/// @param start_time the start time.
/// @param end_time the final time.
/// @param time_delta the initial integration step, it will dynamically change.
/// @return the number of steps it took to perform the integration.
template <class Stepper, class System, class Observer>
constexpr inline auto integrate_adaptive(
    Stepper &stepper,
    Observer &observer,
    System &&system,
    typename Stepper::state_type &state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type time_delta)
{
    using state_type = typename Stepper::state_type;
    // Check if the state vector can (and should) be resized.
    if constexpr (chainsaw::detail::has_resize_v<state_type>) {
        stepper.adjust_size(state);
    }
    // Run until the time reaches the `end_time`, the outer while loop allows to
    // precisely simulate up to end_time. That's why the outer loop is usually
    // simulated 2 times.
    while (chainsaw::detail::less_with_sign(start_time, end_time, time_delta)) {
        // Make sure we don't go beyond the end_time.
        while (chainsaw::detail::less_eq_with_sign(start_time + time_delta, end_time, time_delta)) {
            // Perform one integration step.
            detail::integrate_one_step(stepper, observer, system, state, start_time, time_delta);
            // Advance time.
            start_time += time_delta;
            // Update integration step size.
            time_delta = stepper.get_time_delta();
        }
        // Calculate time step to arrive exactly at end time.
        time_delta = end_time - start_time;
    }
    // Call the observer one last time.
    observer(state, start_time);
    // Return the number of steps it took to integrate.
    return stepper.steps();
}

} // namespace solver