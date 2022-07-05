/// @file stepper_rk4.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

namespace solver
{

template <class State, class Time>
class stepper_rk4 {
public:
    using order_type_t     = unsigned short;
    using time_type_t      = Time;
    using container_type_t = State;
    using value_type_t     = typename State::value_type;

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

    template <class System>
    void do_step(System &system, State &x, const State &dxdt, const Time t, const Time dt)
    {
        const Time val1 = static_cast<Time>(1.0);
        const Time dh   = static_cast<Time>(0.5) * dt;
        const Time th   = t + dh;
        const Time dt6  = dt / static_cast<Time>(6.0);
        const Time dt3  = dt / static_cast<Time>(3.0);
        
        // dt * dxdt = k1 (computed before calling this function)

        // xt = x + dh * dxdt
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, dxdt.begin());

        // dt * m_dxt = k2
        system(m_xt, m_dxt, th);

        // xt = x + dh*m_dxt
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, m_dxt.begin());

        // dt * m_dxm = k3
        system(m_xt, m_dxm, th);

        // xt = x + dt*m_dxm
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dt, m_dxm.begin());

        // dt * m_dxh = k4
        system(m_xt, m_dxh, t + dt);

        // x += dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
        it_algebra::scale_sum_inplace(x.begin(), x.end(), dt6, dxdt.begin(), dt3, m_dxt.begin(), dt3, m_dxm.begin(), dt6, m_dxh.begin());
    }

    template <class System>
    void do_step(System &system, State &x, Time t, Time dt)
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