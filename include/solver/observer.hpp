/// @file observer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-07-05

#pragma once

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