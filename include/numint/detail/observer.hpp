/// @file observer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

#include <cstdint>
#include <iostream>

namespace numint::detail
{

template <class State, class Time>
class Observer {
public:
    inline virtual void operator()(const State &, const Time &)
    {
        // Nothing to do.
    }
};

template <class State, class Time, std::size_t DECIMATION = 1>
class ObserverDecimate : public Observer<State, Time> {
protected:
    explicit ObserverDecimate()
        : decimation_cnt()
    {
        // Nothing to do.
    }

    constexpr bool observe()
    {
        if constexpr (DECIMATION == 0)
            return true;
        if (++decimation_cnt == DECIMATION) {
            decimation_cnt = 0;
            return true;
        }
        return false;
    }

private:
    std::size_t decimation_cnt;
};

template <class State, class Time, std::size_t DECIMATION = 0>
class ObserverPrint : public ObserverDecimate<State, Time, DECIMATION> {
public:
    inline void operator()(const State &x, const Time &t) override
    {
        if (this->observe()) {
            std::cout << t << " " << x << "\n";
        }
    }
};

} // namespace numint::detail
