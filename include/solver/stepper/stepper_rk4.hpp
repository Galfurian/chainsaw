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
        : m_dxdt(),
          m_dxt(),
          m_dxm(),
          m_dxh(),
          m_xt()
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
            m_dxdt.resize(x.size());
            m_dxt.resize(x.size());
            m_dxm.resize(x.size());
            m_dxh.resize(x.size());
            m_xt.resize(x.size());
        }
    }

    template <class System>
    constexpr inline void do_step(System &system, State &x, const State &dxdt, const Time t, const Time dt)
    {
        constexpr Time val1(1.0);
        const Time dh   = static_cast<Time>(0.5) * dt;
        const Time th   = t + dh;
        const Time dt6  = dt / static_cast<Time>(6.0);
        const Time dt3  = dt / static_cast<Time>(3.0);
        // dt * dxdt = k1 (computed before calling this function)
        //system(x, m_dxdt, t);
        // xt = x + dh * dxdt
        detail::it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, dxdt.begin());
        // dt * dxt = k2
        system(m_xt, m_dxt, th);
        // xt = x + dh * dxt
        detail::it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, m_dxt.begin());
        // dt * dxm = k3
        system(m_xt, m_dxm, th);
        // xt = x + dt * dxm
        detail::it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dt, m_dxm.begin());
        // dt * dxh = k4
        system(m_xt, m_dxh, t + dt);
        // x += (dt/6) * (dxdt + 2 * dxt + 2 * dxm + dxh)
        detail::it_algebra::scale_sum_inplace(x.begin(), x.end(), dt6, dxdt.begin(), dt3, m_dxt.begin(), dt3, m_dxm.begin(), dt6, m_dxh.begin());
    }

    template <class System>
    constexpr inline void do_step(System &system, State &x, const Time t, const Time dt)
    {
        // dt * dxdt = k1
        system(x, m_dxdt, t);
        this->do_step(system, x, m_dxdt, t, dt);
    }

private:
    State m_dxdt;
    State m_dxt;
    State m_dxm;
    State m_dxh;
    State m_xt;
};

} // namespace solver