/// @file solver.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "numint/detail/less_with_sign.hpp"
#include "numint/detail/it_algebra.hpp"
#include "numint/detail/type_traits.hpp"

#define NUMINT_MAJOR_VERSION 1 ///< Major version of the library.
#define NUMINT_MINOR_VERSION 0 ///< Minor version of the library.
#define NUMINT_MICRO_VERSION 0 ///< Micro version of the library.

namespace numint
{

namespace detail
{

/// @brief Executes a single integration step and invokes the observer.
///
/// @details This function performs one integration step using the specified
/// stepper and calls the observer with the updated state and time. It is used
/// to advance the state of the system by a single step during numerical
/// integration.
///
/// @tparam Stepper The type of the integration stepper.
/// @tparam System The type of the system being integrated.
/// @tparam Observer The type of the observer function.
///
/// @param stepper The stepper used to perform the integration step.
/// @param observer The observer function to call after the step, receiving the updated state and time.
/// @param system The system being integrated, which defines the equations of motion or dynamics.
/// @param state The current state of the system, which will be updated after the step.
/// @param time The current time at the beginning of the step.
/// @param time_delta The size of the time step for the integration.
template <class Stepper, class System, class Observer>
constexpr inline void integrate_one_step(
    Stepper &stepper,
    Observer &&observer,
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

/// @brief Default termination condition that never ends early.
///
/// @details This function provides a default implementation of a termination
/// condition for numerical integration. It always returns false, indicating
/// that the integration process should not terminate early and should continue
/// until the end time is reached.
///
/// @tparam State The type representing the state of the system.
///
/// @param state The current state of the system (unused).
///
/// @return Always returns false, meaning integration should continue.
template <class State>
constexpr inline bool default_termination_condition(const State &state) noexcept
{
    (void)state;
    return false;
}

} // namespace detail

/// @brief Integrates the system over a fixed time step between the start and end time.
///
/// @details This function performs fixed-step integration of a system over a
/// given time interval using the specified stepper. The observer is invoked
/// after every integration step, and an optional termination condition can be
/// used to stop the integration early.
///
/// @tparam Stepper The type of the integration stepper.
/// @tparam System The type of the system being integrated.
/// @tparam Observer The type of the observer function.
/// @tparam TerminationCondition The type of the termination condition function.
///
/// @param stepper The stepper used to perform the integration.
/// @param observer The observer function to call after each step, receiving the updated state and time.
/// @param system The system being integrated, which defines the equations of motion or dynamics.
/// @param state The initial state of the system, which will be updated during integration.
/// @param start_time The start time for the integration.
/// @param end_time The final time for the integration.
/// @param time_delta The fixed step size for integration.
/// @param check_if_done The termination condition to determine if integration
/// should stop early. Defaults to a function that always returns false.
/// @return The number of steps taken to complete the integration.
template <class Stepper, class System, class Observer, class TerminationCondition = decltype(detail::default_termination_condition<typename Stepper::state_type>)>
constexpr inline auto integrate_fixed(
    Stepper &stepper,
    Observer &&observer,
    System &&system,
    typename Stepper::state_type &state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type time_delta,
    TerminationCondition check_if_done = detail::default_termination_condition<typename Stepper::state_type>) noexcept
{
    using state_type = typename Stepper::state_type;
    // Check if the state vector can (and should) be resized.
    if constexpr (numint::detail::has_resize_v<state_type>) {
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
        // Check if the integration should terminate early by calling the check_if_done function.
        if (check_if_done(state)) {
            break; // Terminate the integration early.
        }
    }
    // Return the number of steps it took to integrate.
    return stepper.steps();
}

/// @brief Integrates the system from the start time to the end time using an
/// adaptive stepper.
///
/// @details This function performs adaptive integration of a system over a
/// specified time interval using the given stepper. It dynamically adjusts the
/// step size to ensure accuracy and stability. The observer is invoked after
/// every integration step, and an optional termination condition can be used to
/// stop the integration early.
///
/// @tparam Stepper The type of the integration stepper.
/// @tparam System The type of the system being integrated.
/// @tparam Observer The type of the observer function.
/// @tparam TerminationCondition The type of the termination condition function.
///
/// @param stepper The stepper used to perform the integration.
/// @param observer The observer function to call after each step, receiving the updated state and time.
/// @param system The system being integrated, which defines the equations of motion or dynamics.
/// @param state The initial state of the system, which will be updated during integration.
/// @param start_time The start time for the integration.
/// @param end_time The final time for the integration.
/// @param time_delta The initial step size for integration. This may be dynamically adjusted.
/// @param check_if_done The termination condition to determine if integration
/// should stop early. Defaults to a function that always returns false.
///
/// @return The number of steps taken to complete the integration.
template <class Stepper, class System, class Observer, class TerminationCondition = decltype(detail::default_termination_condition<typename Stepper::state_type>)>
constexpr inline auto integrate_adaptive(
    Stepper &stepper,
    Observer &&observer,
    System &&system,
    typename Stepper::state_type &state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type time_delta,
    TerminationCondition check_if_done = detail::default_termination_condition<typename Stepper::state_type>)
{
    using state_type = typename Stepper::state_type;
    
    // Adjust the stepper's internal size if the state supports resizing.
    if constexpr (numint::detail::has_resize_v<state_type>) {
        stepper.adjust_size(state);
    }

    // Run until the time reaches the `end_time`, the outer while loop allows to
    // precisely simulate up to end_time. That's why the outer loop is usually
    // simulated 2 times.
    while (numint::detail::less_with_sign(start_time, end_time, time_delta)) {
        // Make sure we don't go beyond the end_time.
        while (numint::detail::less_eq_with_sign(start_time + time_delta, end_time, time_delta)) {
            // Perform one integration step.
            detail::integrate_one_step(stepper, observer, system, state, start_time, time_delta);
            // Advance time.
            start_time += time_delta;
            // Update integration step size.
            time_delta = stepper.get_time_delta();
            // Check if the integration should terminate early by calling the check_if_done function.
            if (check_if_done(state)) {
                break; // Terminate the integration early.
            }
        }
        // Calculate time step to arrive exactly at end time.
        time_delta = end_time - start_time;
    }
    // Call the observer one last time.
    observer(state, start_time);
    // Return the number of steps it took to integrate.
    return stepper.steps();
}

} // namespace numint