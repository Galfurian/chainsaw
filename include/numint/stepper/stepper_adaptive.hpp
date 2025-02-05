/// @file stepper_adaptive.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "numint/detail/it_algebra.hpp"
#include "numint/detail/type_traits.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>

namespace numint
{

/// Types of truncation error formulas.
enum class ErrorFormula : unsigned char {
    Absolute, ///< Use the absolute truncation error.
    Relative, ///< Use the relative truncation error.
    Mixed     ///< Use a mixed absolute and relative truncation error.
};

/// @brief It dynamically controlls the step-size of a stepper.
/// @tparam Stepper The stepper we rely upon.
/// @tparam Iterations The number of iterations we are going to do while
/// integrating, higher values means more accurate results, but computationally
/// expensive.
/// @tparam Error The type of error formula we rely upon.
template <class Stepper, int Iterations = 2, ErrorFormula Error = ErrorFormula::Absolute>
class stepper_adaptive
{
public:
    /// @brief Type of internal fixed-step stepper we are using.
    using stepper_type                        = Stepper;
    /// @brief Type used for the order of the stepper.
    using order_type                          = typename Stepper::order_type;
    /// @brief Type used to keep track of time.
    using time_type                           = typename Stepper::time_type;
    /// @brief The state vector.
    using state_type                          = typename Stepper::state_type;
    /// @brief Type of value contained in the state vector.
    using value_type                          = typename Stepper::state_type::value_type;
    /// @brief Determines if this is an adaptive stepper or not.
    static constexpr bool is_adaptive_stepper = true;

    /// @brief Creates a new adaptive stepper.
    stepper_adaptive()
        : m_stepper_main()
        , m_stepper_tuner()
        , m_tollerance(0.0001)
        , m_time_delta(1e-12)
        , m_min_delta(1e-12)
        , m_max_delta(1)
        , m_t_err(.0)
        , m_t_err_abs(.0)
        , m_t_err_rel(.0)

    {
        // Nothing to do.
    }

    /// @brief Destructor.
    ~stepper_adaptive() = default;

    /// @brief Copy construction is disabled.
    stepper_adaptive(const stepper_adaptive &other) = delete;

    /// @brief Copy assignment is disabled.
    auto operator=(const stepper_adaptive &other) -> stepper_adaptive & = delete;

    /// @brief Move constructor.
    stepper_adaptive(stepper_adaptive &&other) noexcept = default;

    /// @brief Move assignment operator.
    auto operator=(stepper_adaptive &&other) noexcept -> stepper_adaptive & = default;

    /// @brief Sets the tolerance for step-size control.
    ///
    /// @param tollerance The tolerance value to use for adjusting the step size.
    constexpr void set_tollerance(value_type tollerance) { m_tollerance = tollerance; }

    /// @brief Sets the minimum allowed step size.
    ///
    /// @param min_delta The minimum step size.
    constexpr void set_min_delta(value_type min_delta) { m_min_delta = min_delta; }

    /// @brief Sets the maximum allowed step size.
    ///
    /// @param max_delta The maximum step size.
    constexpr void set_max_delta(value_type max_delta) { m_max_delta = max_delta; }

    /// @brief The order of the stepper we rely upon.
    /// @return the order of the internal stepper.
    constexpr auto order_step() const -> order_type { return m_stepper_main.order_step(); }

    /// @brief Retrieves the current adaptive step size.
    /// @return The current step size as a `time_type` value.
    constexpr auto get_time_delta() const -> time_type { return m_time_delta; }

    /// @brief Adjusts the size of the internal state vectors.
    /// @param reference a reference state vector vector.
    void adjust_size(const state_type &reference)
    {
        m_stepper_main.adjust_size(reference);
        m_stepper_tuner.adjust_size(reference);
    }

    /// @brief Returns the number of steps the stepper executed up until now.
    /// @return the number of integration steps.
    constexpr auto steps() const { return m_steps; }

    /// @brief Performs one integration step using the provided system.
    ///
    /// @details This function advances the state of the system by one step
    /// using the integration method described below:
    ///
    /// Compute values of (0):
    ///     y_{n + 1} = y_n + h * f(t_n, y_n)
    ///
    /// Compute values of (1):
    ///     y_{n + 0.5} = y_n + 0.5 * h * f(t_n, y_n)
    ///     y_{n + 1} = y_{n + 0.5} + 0.5 * h * f(t_n, y_n)
    ///
    /// @tparam System The type of the system being integrated.
    ///
    /// @param system The system that defines the equations of motion or dynamics.
    /// @param x The state of the system, which will be updated after this step.
    /// @param t The current time.
    /// @param dt The time step to use for the integration.
    template <class System>
    constexpr void do_step(System &&system, state_type &x, const time_type t, const time_type dt)
    {
        using detail::it_algebra::max_abs_diff;
        using detail::it_algebra::max_comb_diff;
        using detail::it_algebra::max_rel_diff;

        // Copy the step size.
        m_time_delta = dt;
        // Copy the initial state.
        state_type y(x);
        // Compute values of (0).
        m_stepper_main.do_step(std::forward<System>(system), y, t, m_time_delta);
        // Compute values of (1).
        if constexpr (Iterations <= 2) {
            const time_type dh = m_time_delta * .5;
            m_stepper_tuner.do_step(std::forward<System>(system), x, t + dh, dh);
            m_stepper_tuner.do_step(std::forward<System>(system), x, t + m_time_delta, dh);
        } else {
            const time_type dh = m_time_delta * (1. / Iterations);
            for (unsigned i = 0; i < Iterations; ++i) {
                m_stepper_tuner.do_step(std::forward<System>(system), x, t + (dh * i), dh);
            }
        }
        // Calculate truncation error.
        if constexpr (Error == ErrorFormula::Absolute) {
            // Get absolute truncation error.
            m_t_err_abs = max_abs_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            m_time_delta *= 0.9 * std::min(std::max(std::pow(m_tollerance / (2 * m_t_err_abs), 0.2), 0.3), 2.);
        } else if constexpr (Error == ErrorFormula::Relative) {
            // Get relative truncation error.
            m_t_err_rel = max_rel_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            m_time_delta *= 0.9 * std::min(std::max(std::pow(m_tollerance / (2 * m_t_err_rel), 0.2), 0.3), 2.);
        } else {
            // Get mixed truncation error.
            m_t_err = max_comb_diff<value_type>(x.begin(), x.end(), y.begin(), y.end());
            // Update the time-delta.
            m_time_delta *= 0.9 * std::min(std::max(std::pow(m_tollerance / (2 * m_t_err), 0.2), 0.3), 2.);
        }
        // Check boundaries.
        m_time_delta = std::min(std::max(m_time_delta, m_min_delta), m_max_delta);
        // Increase the number of steps.
        ++m_steps;

        /*
        // Copy the step size.
        m_time_delta = dt;
        // Copy the initial state.
        state_type y0, y1;
        while (true) {
            y0 = x, y1 = x;

            // Compute values of (0).
            m_stepper_main.do_step(std::forward<System>(system), y0, t, m_time_delta);

            // Compute values of (1).
            if constexpr (Iterations <= 2) {
                const time_type dh = m_time_delta * .5;
                m_stepper_tuner.do_step(std::forward<System>(system), y1, t + dh, dh);
                m_stepper_tuner.do_step(std::forward<System>(system), y1, t + m_time_delta, dh);
            } else {
                const time_type dh = m_time_delta * (1. / Iterations);
                for (unsigned i = 0; i < Iterations; ++i) {
                    m_stepper_tuner.do_step(std::forward<System>(system), y1, t + dh * i, dh);
                }
            }

            // Calculate truncation error.
            if constexpr (Error == ErrorFormula::Absolute) {
                // Get absolute truncation error.
                m_t_err = max_abs_diff<value_type>(y0.begin(), y0.end(), y1.begin(), y1.end());
            } else if constexpr (Error == ErrorFormula::Relative) {
                // Get relative truncation error.
                m_t_err = max_rel_diff<value_type>(y0.begin(), y0.end(), y1.begin(), y1.end());
            } else {
                // Get mixed truncation error.
                m_t_err = max_comb_diff<value_type>(y0.begin(), y0.end(), y1.begin(), y1.end());
            }

            // Update the time-delta.
            m_time_delta *= 0.9 * std::min(std::max(std::pow(m_tollerance / (2 * m_t_err), 0.2), 0.3), 2.);
            // Check boundaries.
            m_time_delta = std::min(std::max(m_time_delta, m_min_delta), m_max_delta);

            // Increase the number of steps.
            ++m_steps;

            // Check if error is within tolerance
            if (m_t_err <= m_tollerance) {
                std::cout << m_t_err << "\n";
                // Update the state.
                x = y1;
                break;
            }
        }
        */
    }

private:
    /// The main stepper.
    stepper_type m_stepper_main;
    /// A temporary stepper we use to tune the main stepper.
    stepper_type m_stepper_tuner;
    /// The tollerance value we use to tune the step-size.
    time_type m_tollerance;
    /// A copy of the step-size.
    time_type m_time_delta;
    /// The minimum step-size.
    time_type m_min_delta;
    /// The maximum step-size.
    time_type m_max_delta;
    /// Holds the error between the main stepper and the temporary stepper.
    value_type m_t_err;
    /// Holds the absolute error between the main stepper and the temporary stepper.
    value_type m_t_err_abs;
    /// Holds the relative error between the main stepper and the temporary stepper.
    value_type m_t_err_rel;
    /// The number of steps of integration.
    uint64_t m_steps{};
};

} // namespace numint
