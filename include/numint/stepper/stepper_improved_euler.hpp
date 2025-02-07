/// @file stepper_improved_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "numint/detail/it_algebra.hpp"
#include "numint/detail/type_traits.hpp"

namespace numint
{

/// @brief Stepper implementing Heun's method for numerical integration (also
/// known as the Improved Euler Method).
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_improved_euler
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
    stepper_improved_euler() = default;

    /// @brief Destructor.
    ~stepper_improved_euler() = default;

    /// @brief Copy constructor.
    /// @param other The logger instance to copy from.
    stepper_improved_euler(const stepper_improved_euler &other) = delete;

    /// @brief Move constructor.
    /// @param other The logger instance to move from.
    stepper_improved_euler(stepper_improved_euler &&other) noexcept = default;

    /// @brief Copy assignment operator.
    /// @param other The logger instance to copy from.
    /// @return Reference to the logger instance.
    auto operator=(const stepper_improved_euler &other) -> stepper_improved_euler & = delete;

    /// @brief Move assignment operator.
    /// @param other The logger instance to move from.
    /// @return Reference to the logger instance.
    auto operator=(stepper_improved_euler &&other) noexcept -> stepper_improved_euler & = default;

    /// @brief Returns the order of the stepper.
    /// @return The order of the internal stepper, which is 2 for the Improved Euler method.
    constexpr auto order_step() const -> order_type { return 1; }

    /// @brief Adjusts the size of the internal state vectors based on a reference.
    /// @param reference A reference state vector used for size adjustment.
    void adjust_size(const state_type &reference)
    {
        if constexpr (numint::detail::has_resize<state_type>::value) {
            m_dxdt1.resize(reference.size()); // Resize m_dxdt1 if supported.
            m_dxdt2.resize(reference.size()); // Resize m_dxdt2 if supported.
            m_x.resize(reference.size());     // Resize m_x if supported.
        }
    }

    /// @brief Returns the number of steps executed by the stepper so far.
    /// @return The number of integration steps executed.
    constexpr auto steps() const { return m_steps; }

    /// @brief Performs a single integration step using Heun's method (Improved Euler method).
    /// @param system The system to integrate.
    /// @param x The initial state vector.
    /// @param t The initial time.
    /// @param dt The time step for integration.
    template <class System>
    constexpr void do_step(System &&system, state_type &x, const time_type t, const time_type dt) noexcept
    {
        // Calculate the derivative at the initial point:
        //      dxdt1 = system(x, t);
        std::forward<System>(system)(x, m_dxdt1, t);

        // Calculate the state at the next time point using Euler's method:
        //      m_x(t + dt) = x(t) + dxdt1 * dt;
        detail::it_algebra::sum_operation(
            m_x.begin(), m_x.end(), std::multiplies<>(), 1., x.begin(), dt, m_dxdt1.begin());

        // Calculate the derivative at the midpoint:
        //      dxdt2 = system(m_x, t + dt);
        std::forward<System>(system)(m_x, m_dxdt2, t + dt);

        // Update the state vector using the average of the derivatives:
        //      x(t + dt) = x(t) + (dt / 2) * (dxdt1 + dxdt2);
        detail::it_algebra::accumulate_operation(
            x.begin(), x.end(), std::multiplies<>(), dt * .5, m_dxdt1.begin(), dt * .5, m_dxdt2.begin());

        // Increment the number of integration steps.
        ++m_steps;
    }

private:
    /// Keeps track of the first derivative of the state.
    state_type m_dxdt1;

    /// Keeps track of the second derivative of the state.
    state_type m_dxdt2;

    /// Temporary state vector for intermediate calculations.
    state_type m_x;

    /// The number of steps taken during integration.
    unsigned long m_steps{};
};

} // namespace numint
