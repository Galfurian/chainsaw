/// @file stepper_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "solver/detail/type_traits.hpp"
#include "solver/detail/it_algebra.hpp"

namespace solver
{

template <class State, class Time>
class stepper_euler {
public:
    using order_type_t = unsigned short;
    using time_type_t  = Time;
    using state_type_t = State;
    using value_type_t = typename State::value_type;

    stepper_euler()
        : m_dxdt(),
          m_steps()
    {
        // Nothing to do.
    }

    stepper_euler(const stepper_euler &other) = delete;

    stepper_euler &operator=(const stepper_euler &other) = delete;

    constexpr inline order_type_t order_step() const
    {
        return 1;
    }

    void adjust_size(state_type_t &x)
    {
        if constexpr (solver::detail::has_resize<state_type_t>::value) {
            m_dxdt.resize(x.size());
        }
    }

    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Performs one step with the knowledge of dxdt(t)
    /// @param system
    /// @param x
    /// @param dxdt
    /// @param t
    /// @param dt
    template <class System>
    constexpr void do_step(System &, State &x, const State &dxdt, Time, Time dt) noexcept
    {
        // x(t + 1) = x(t) + dxdt * dt;
        detail::it_algebra::increment(x.begin(), x.end(), dxdt.begin(), dt);
        // Increase the number of steps.
        ++m_steps;
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        // dxdt = system(x, t);
        system(x, m_dxdt, t);
        // x(t + 1) = x(t) + dxdt * dt;
        this->do_step(system, x, m_dxdt, t, dt);
    }

private:
    State m_dxdt;
    unsigned long m_steps;
};

} // namespace solver