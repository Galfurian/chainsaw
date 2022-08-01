/// @file stepper_rk4.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "solver/detail/type_traits.hpp"
#include "solver/detail/it_algebra.hpp"

namespace solver
{

template <class State, class Time>
class stepper_rk4 {
public:
    using order_type_t = unsigned short;
    using time_type_t  = Time;
    using state_type_t = State;
    using value_type_t = typename State::value_type;

    stepper_rk4()
        : m_dxdt1(),
          m_dxdt2(),
          m_dxdt3(),
          m_dxdt4(),
          m_x()
    {
        // Nothing to do.
    }

    constexpr inline order_type_t order_step() const
    {
        return 4;
    }

    constexpr inline void adjust_size(state_type_t &x)
    {
        if constexpr (solver::detail::has_resize<state_type_t>::value) {
            m_dxdt1.resize(x.size());
            m_dxdt2.resize(x.size());
            m_dxdt3.resize(x.size());
            m_dxdt4.resize(x.size());
            m_x.resize(x.size());
        }
    }

    template <class System>
    constexpr inline void do_step(System &system, State &x, const State &dxdt, const Time t, const Time dt)
    {
        const Time dt6 = dt / static_cast<Time>(6);
        const Time dt3 = dt / static_cast<Time>(3);
        const Time dh  = dt / static_cast<Time>(2);
        const Time th  = t + dh;

        // dxdt    = f(x, t) (computed before calling this function)
        // m_x += x + dh * dxdt
        detail::it_algebra::scale_sum(
            m_x.begin(), m_x.end(),
            x.begin(),
            dh, dxdt.begin());

        // m_dxdt2 = f(m_x, t + dh)
        // m_x += x + dh * m_dxdt2
        system(m_x, m_dxdt2, th);
        detail::it_algebra::scale_sum(
            m_x.begin(), m_x.end(),
            x.begin(),
            dh, m_dxdt2.begin());

        // m_dxdt3 = f(m_x, t + dh)
        // m_x += x + dt * m_dxdt3
        system(m_x, m_dxdt3, th);
        detail::it_algebra::scale_sum(
            m_x.begin(), m_x.end(),
            x.begin(),
            dt, m_dxdt3.begin());

        // m_dxdt4 = f(m_x, t + dt)
        // x += (dt/6)*dxdt + (dt/3)*m_dxdt2 + (dt/3)*m_dxdt3 + (dt/6)*m_dxdt4
        system(m_x, m_dxdt4, t + dt);
        detail::it_algebra::scale_sum_inplace(
            x.begin(), x.end(),
            dt6, dxdt.begin(),
            dt3, m_dxdt2.begin(),
            dt3, m_dxdt3.begin(),
            dt6, m_dxdt4.begin());
    }

    template <class System>
    constexpr inline void do_step(System &system, State &x, const Time t, const Time dt)
    {
        // dt * dxdt = f(x, t)
        system(x, m_dxdt1, t);
        this->do_step(system, x, m_dxdt1, t, dt);
    }

private:
    State m_dxdt1;
    State m_dxdt2;
    State m_dxdt3;
    State m_dxdt4;
    State m_x;
};

} // namespace solver