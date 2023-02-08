/// @file observer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#pragma once

namespace solver::detail
{

template <std::size_t DECIMATION>
class DecimationObserver {
protected:
    explicit DecimationObserver()
        : decimation_cnt()
    {
        // Nothing to do.
    }

    constexpr inline bool observe()
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

class NoObserver {
public:
    template <class State, class Time>
    inline void operator()(const State &, const Time &)
    {
        // Nothing to do.
    }
};

template <std::size_t DECIMATION = 0>
class ObserverPrint : public DecimationObserver<DECIMATION> {
public:
    template <class State, class Time>
    inline void operator()(const State &x, const Time &t)
    {
        if (this->observe())
            std::cout << t << " " << x << "\n";
    }
};

} // namespace solver::detail
