/// @file stepper_euler.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "chainsaw/detail/type_traits.hpp"
#include "chainsaw/detail/it_algebra.hpp"

namespace chainsaw
{

/// @brief Stepper implementing the Euler method for numerical integration.
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class stepper_euler {
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
        : m_dxdt(), ///< Initializes the derivative state vector.
          m_steps() ///< Initializes the step count.
    {
        // Nothing to do.
    }

    /// @brief Deleted copy constructor.
    stepper_euler(const stepper_euler &other) = delete;

    /// @brief Deleted copy assignment operator.
    stepper_euler &operator=(const stepper_euler &other) = delete;

    /// @brief Returns the order of the stepper.
    /// @return The order of the internal stepper, which is 1 for Euler method.
    constexpr inline order_type order_step() const
    {
        return 1;
    }

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
    constexpr inline auto steps() const
    {
        return m_steps;
    }

    /// @brief Performs a single integration step using Euler's method.
    /// @tparam System The type of the system representing the differential equations.
    /// @param system The system to integrate.
    /// @param x The initial state vector.
    /// @param t The initial time.
    /// @param dt The time step for integration.
    template <class System>
    constexpr void do_step(System &&system, state_type &x, const time_type t, const time_type dt) noexcept
    {
        // Calculate the derivative at the current time.
        system(x, m_dxdt, t);

        // Update the state vector using Euler's method:
        //      x(t + dt) = x(t) + dxdt * dt.
        detail::it_algebra::scale_accumulate(x.begin(), x.end(), m_dxdt.begin(), dt);

        // Increment the number of integration steps.
        ++m_steps;
    }

private:
    /// Keeps track of the derivative of the state.
    state_type m_dxdt;

    /// The number of steps taken during integration.
    unsigned long m_steps;
};

} // namespace chainsaw
