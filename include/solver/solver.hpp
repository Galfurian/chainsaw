/// @file solver.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "detail/less_with_sign.hpp"
#include "detail/it_algebra.hpp"

namespace solver
{

namespace detail
{

template <class Stepper, class System, class Observer>
constexpr inline void integrate_one_step_const(
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
constexpr inline auto integrate_const(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta) noexcept
{
    if constexpr (solver::detail::has_resize<typename Stepper::state_type_t>::value) {
        stepper.adjust_size(state);
    }
    std::size_t iteration = 0;
    observer(state, start_time);
    while (start_time < end_time) {
        detail::integrate_one_step_const(stepper, observer, system, state, start_time, time_delta);
        start_time += time_delta;
        ++iteration;
    }
    return iteration;
}

template <class Stepper, class System, class Observer>
constexpr inline auto integrate_adaptive_simple(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta) noexcept
{
    using state_type_t = typename Stepper::state_type_t;

    // Check if the state vector can (and should) be resized.
    if constexpr (solver::detail::has_resize<state_type_t>::value) {
        stepper.adjust_size(state);
    }

    // Integrate with 
    std::size_t iteration = integrate_const(stepper, observer, system, state, start_time, end_time, time_delta);
    // Make a last step to end exactly at end_time.
    if (less_with_sign(start_time + time_delta * iteration, end_time, time_delta)) {
        detail::integrate_one_step_const(stepper, observer, system, state, start_time, time_delta);
        ++iteration;
    }
    return iteration;
}

template <class Stepper, class System, class Observer>
constexpr inline auto integrate_adaptive(
    Stepper &stepper,
    Observer &observer,
    System &system,
    typename Stepper::state_type_t &state,
    typename Stepper::time_type_t start_time,
    typename Stepper::time_type_t end_time,
    typename Stepper::time_type_t time_delta)
{
    using state_type_t = typename Stepper::state_type_t;
    using time_type_t = typename Stepper::time_type_t;
    // Check if the state vector can (and should) be resized.
    if constexpr (solver::detail::has_resize<state_type_t>::value) {
        stepper.adjust_size(state);
    }

    // Initialize the number of iterations.
    std::size_t iteration = 0;

    // Initilize the stepper.
    stepper.initialize(state, start_time, time_delta);

    // Run until the time reaches the `end_time`.
    while (solver::detail::less_with_sign(stepper.current_time(), end_time, stepper.current_time_step())) {
        std::cout << stepper.current_time() << "\n";
        // Make sure we don't go beyond the end_time.
        while (solver::detail::less_eq_with_sign(static_cast<time_type_t>(stepper.current_time() + stepper.current_time_step()), end_time, stepper.current_time_step())) {
            stepper.do_step(system);
            observer(stepper.current_state(), stepper.current_time());
            ++iteration;
        }
        // Calculate time step to arrive exactly at end time.
        stepper.initialize(stepper.current_state(), stepper.current_time(), static_cast<time_type_t>(end_time - stepper.current_time()));
    }
    observer(stepper.current_state(), stepper.current_time());
    // Overwrite state with the final point.
    state = stepper.current_state();
    return iteration;
}

} // namespace solver