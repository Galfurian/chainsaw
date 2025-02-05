/// @file observer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include <cstdint>
#include <iostream>

namespace numint::detail
{

/// @brief Observer class.
///
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
template <class State, class Time>
class Observer
{
public:
    /// @brief Default constructor.
    Observer() = default;

    /// @brief Copy constructor.
    Observer(const Observer &other) = default;

    /// @brief Copy assignment operator.
    auto operator=(const Observer &other) -> Observer & = default;

    /// @brief Move constructor.
    Observer(Observer &&other) noexcept = default;

    /// @brief Move assignment operator.
    auto operator=(Observer &&other) noexcept -> Observer & = default;

    /// @brief Destructor.
    virtual ~Observer() = default;

    /// @brief Perform the observation.
    /// @param x The state vector.
    /// @param t The time.
    virtual void operator()(const State &x, const Time &t) { (void)x, (void)t; }
};

/// @brief Observer class that decimates the observation.
///
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
/// @tparam DECIMATION The decimation factor.
template <class State, class Time, std::size_t DECIMATION = 1>
class ObserverDecimate : public Observer<State, Time>
{
protected:
    /// @brief Constructor.
    explicit ObserverDecimate()

    {
        // Nothing to do.
    }

    /// @brief Determines if the observer should observe the current state.
    /// @return true if the observer should observe the current state, false otherwise.
    constexpr auto observe() -> bool
    {
        if constexpr (DECIMATION == 0) {
            return true;
        }
        if (++decimation_cnt == DECIMATION) {
            decimation_cnt = 0;
            return true;
        }
        return false;
    }

private:
    /// @brief The decimation counter.
    std::size_t decimation_cnt{};
};

/// @brief Observer that prints the state vector.
/// @tparam State The state vector type.
/// @tparam Time The datatype used to hold time.
/// @tparam DECIMATION The decimation factor.
template <class State, class Time, std::size_t DECIMATION = 0>
class ObserverPrint : public ObserverDecimate<State, Time, DECIMATION>
{
public:
    void operator()(const State &x, const Time &t) override
    {
        if (this->observe()) {
            std::cout << t << " " << x << "\n";
        }
    }
};

} // namespace numint::detail
