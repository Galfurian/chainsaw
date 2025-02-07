/// @file stepper_rk4.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "numint/detail/it_algebra.hpp"
#include "numint/detail/type_traits.hpp"

namespace numint
{

/// @brief Stepper implementing a fourth-order Runge-Kutta method.
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_rk4
{
public:
    /// @brief Type used for the order of the stepper.
    using order_type = unsigned short;

    /// @brief Type used to keep track of time.
    using time_type = Time;

    /// @brief The state vector type.
    using state_type = State;

    /// @brief Type of value contained in the state vector.
    using value_type = typename state_type::value_type;

    /// @brief Indicates whether this is an adaptive stepper.
    static constexpr bool is_adaptive_stepper = false;

    /// @brief Constructs a new stepper.
    stepper_rk4() = default;

    /// @brief Destructor.
    ~stepper_rk4() = default;

    /// @brief Copy constructor.
    /// @param other The logger instance to copy from.
    stepper_rk4(const stepper_rk4 &other) = delete;

    /// @brief Move constructor.
    /// @param other The logger instance to move from.
    stepper_rk4(stepper_rk4 &&other) noexcept = default;

    /// @brief Copy assignment operator.
    /// @param other The logger instance to copy from.
    /// @return Reference to the logger instance.
    auto operator=(const stepper_rk4 &other) -> stepper_rk4 & = delete;

    /// @brief Move assignment operator.
    /// @param other The logger instance to move from.
    /// @return Reference to the logger instance.
    auto operator=(stepper_rk4 &&other) noexcept -> stepper_rk4 & = default;

    /// @brief Returns the order of the stepper.
    /// @return The order of the internal stepper, which is 4 for RK4.
    constexpr auto order_step() const -> order_type { return 4; }

    /// @brief Adjusts the size of the internal state vectors based on a reference.
    /// @param reference A reference state vector used for size adjustment.
    constexpr void adjust_size(const state_type &reference)
    {
        if constexpr (detail::has_resize<state_type>::value) {
            m_dxdt1.resize(reference.size());
            m_dxdt2.resize(reference.size());
            m_dxdt3.resize(reference.size());
            m_dxdt4.resize(reference.size());
            m_x.resize(reference.size());
        }
    }

    /// @brief Returns the number of steps executed by the stepper so far.
    /// @return The number of integration steps executed.
    constexpr auto steps() const { return m_steps; }

    /// @brief Performs a single integration step using the fourth-order Runge-Kutta method.
    /// @tparam System The type of the system representing the differential equations.
    /// @param system The system to integrate.
    /// @param x The initial state vector.
    /// @param t The initial time.
    /// @param dt The time step for integration.
    template <class System>
    constexpr void do_step(System &&system, state_type &x, const time_type t, const time_type dt) noexcept
    {
        // Here is the idea:
        //  - m_dxdt1 : Slope at the beginning of the interval
        //  - m_dxdt2 : Slope at the midpoint of the interval
        //  - m_dxdt3 : Another slope at the midpoint of the interval
        //  - m_dxdt4 : Slope at the end of the interval

        // Step 1: Calculate the slope at the beginning of the interval (m_dxdt1):
        //      m_dxdt1 = f(x, t);
        std::forward<System>(system)(x, m_dxdt1, t);

        // Update temporary state using the slope at the beginning and move halfway forward:
        //      m_x(t + dt * 0.5) = x(t) + m_dxdt1 * dt * 0.5;
        detail::it_algebra::sum_operation(
            m_x.begin(), m_x.end(), std::multiplies<>(), 1.0, x.begin(), 0.5 * dt, m_dxdt1.begin());

        // Step 2: Calculate the slope at the midpoint of the interval (m_dxdt2):
        //      m_dxdt2 = f(m_x, t + 0.5 * dt);
        std::forward<System>(system)(m_x, m_dxdt2, t + (0.5 * dt));

        // Update temporary state using the slope at the midpoint and move halfway forward again:
        //      m_x(t + dt * 0.5) = x(t) + m_dxdt2 * dt * 0.5;
        detail::it_algebra::sum_operation(
            m_x.begin(), m_x.end(), std::multiplies<>(), 1.0, x.begin(), 0.5 * dt, m_dxdt2.begin());

        // Step 3: Calculate another slope at the midpoint of the interval (m_dxdt3):
        //      m_dxdt3 = f(m_x, t + 0.5 * dt);
        std::forward<System>(system)(m_x, m_dxdt3, t + (0.5 * dt));

        // Update temporary state using the slope at the midpoint and move to the end of the interval:
        //      m_x(t + dt) = x(t) + m_dxdt3 * dt;
        detail::it_algebra::sum_operation(
            m_x.begin(), m_x.end(), std::multiplies<>(), 1.0, x.begin(), dt, m_dxdt3.begin());

        // Step 4: Calculate the slope at the end of the interval (m_dxdt4):
        //      m_dxdt4 = f(m_x, t + dt);
        std::forward<System>(system)(m_x, m_dxdt4, t + dt);

        // Update each component of the state vector using the weighted average
        // of the slopes: m_dxdt1, m_dxdt2, m_dxdt3, and m_dxdt4.
        detail::it_algebra::accumulate_operation(
            x.begin(), x.end(), std::multiplies<>(), dt * (1. / 6.), m_dxdt1.begin(), dt * (2. / 6.), m_dxdt2.begin(),
            dt * (2. / 6.), m_dxdt3.begin(), dt * (1. / 6.), m_dxdt4.begin());

        // Increase the number of steps.
        ++m_steps;
    }

private:
    /// Support vectors for the slopes.
    state_type m_dxdt1, m_dxdt2, m_dxdt3, m_dxdt4, m_x;

    /// The number of steps of integration.
    unsigned long m_steps{};
};

} // namespace numint
