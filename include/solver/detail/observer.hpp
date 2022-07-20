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

} // namespace solver::detail
