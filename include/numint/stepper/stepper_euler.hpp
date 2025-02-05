/// @file stepper_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "numint/detail/it_algebra.hpp"
#include "numint/detail/type_traits.hpp"

namespace numint
{

/// @brief Stepper implementing the Euler method for numerical integration.
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_euler
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
    stepper_euler()
        : m_dxdt()
        
    {
        // Nothing to do.
    }

    /// @brief Destructor.
    ~stepper_euler() = default;

    /// @brief Copy construction is disabled.
    stepper_euler(const stepper_euler &other) = delete;

    /// @brief Copy assignment is disabled.
    auto operator=(const stepper_euler &other) -> stepper_euler & = delete;

    /// @brief Move constructor.
    stepper_euler(stepper_euler &&other) noexcept = default;

    /// @brief Move assignment operator.
    auto operator=(stepper_euler &&other) noexcept -> stepper_euler & = default;

    /// @brief Returns the order of the stepper.
    /// @return The order of the internal stepper, which is 1 for Euler method.
    constexpr auto order_step() const -> order_type { return 1; }

    /// @brief Adjusts the size of the internal state vector based on a reference.
    /// @param reference A reference state vector used for size adjustment.
    void adjust_size(const state_type &reference)
    {
        if constexpr (detail::has_resize<state_type>::value) {
            m_dxdt.resize(reference.size()); // Resize m_dxdt if supported.
        }
    }

    /// @brief Returns the number of steps executed by the stepper so far.
    /// @return The number of integration steps executed.
    constexpr auto steps() const { return m_steps; }

    /// @brief Performs a single integration step using Euler's method.
    /// @tparam System The type of the system representing the differential equations.
    /// @param system The system to integrate.
    /// @param x The initial state vector.
    /// @param t The initial time.
    /// @param dt The time step for integration.
    template <class System>
    void do_step(System &&system, state_type &x, const time_type t, const time_type dt)
    {
        // Calculate the derivative at the current time.
        std::forward<System>(system)(x, m_dxdt, t);

        // Update the state vector using Euler's method:
        //      x(t + dt) = x(t) + dxdt * dt.
        detail::it_algebra::accumulate_operation(x.begin(), x.end(), std::multiplies<>(), dt, m_dxdt.begin());

        // Increment the number of integration steps.
        ++m_steps;
    }

private:
    /// Keeps track of the derivative of the state.
    state_type m_dxdt;

    /// The number of steps taken during integration.
    unsigned long m_steps{};
};

} // namespace numint
