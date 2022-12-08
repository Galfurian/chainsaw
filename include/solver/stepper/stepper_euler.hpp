/// @file stepper_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "solver/detail/type_traits.hpp"
#include "solver/detail/it_algebra.hpp"

namespace solver
{

/// @brief Stepper implementing Euler method.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_euler {
public:
    using order_type_t = unsigned short;
    using time_type_t  = Time;
    using state_type_t = State;
    using value_type_t = typename State::value_type;

    /// @brief Creates a new stepper.
    stepper_euler()
        : m_dxdt(),
          m_steps()
    {
        // Nothing to do.
    }

    /// @brief Nope.
    stepper_euler(const stepper_euler &other) = delete;

    /// @brief Nope.
    stepper_euler &operator=(const stepper_euler &other) = delete;

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr inline order_type_t order_step() const
    {
        return 1;
    }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type_t &reference)
    {
        if constexpr (solver::detail::has_resize<state_type_t>::value) {
            m_dxdt.resize(reference.size());
        }
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Integrates on step.
    /// @param system the system we are integrating.
    /// @param x the initial state.
    /// @param t the initial time.
    /// @param dt the step-size.
    template <class System>
    constexpr void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        // dxdt = system(x, t);
        system(x, m_dxdt, t);
        // x(t + 1) = x(t) + dxdt * dt;
        detail::it_algebra::increment(x.begin(), x.end(), m_dxdt.begin(), dt);
        // Increase the number of steps.
        ++m_steps;
    }

private:
    /// Keeps track of state evolution.
    State m_dxdt;
    /// The number of steps of integration.
    unsigned long m_steps;
};

} // namespace solver