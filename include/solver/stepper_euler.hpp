/// @file stepper_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include "type_traits.hpp"

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
        : m_dxdt()
    {
        // Nothing to do.
    }

    constexpr inline order_type_t order_step() const
    {
        return 1;
    }

    void adjust_size(state_type_t &x)
    {
        if constexpr (has_resize<state_type_t>::value) {
            m_dxdt.resize(x.size());
        }
    }

    /// @brief Performs one step with the knowledge of dxdt(t)
    /// @param system
    /// @param x
    /// @param dxdt
    /// @param t
    /// @param dt
    template <class System>
    constexpr void do_step(System &system, State &x, const State &dxdt, Time t, Time dt) noexcept
    {
        it_algebra::increment(x.begin(), x.end(), dxdt.begin(), dt);
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        system(x, m_dxdt, t);
        this->do_step(system, x, m_dxdt, t, dt);
    }

private:
    State m_dxdt;
};

} // namespace solver