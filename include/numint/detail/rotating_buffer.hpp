/// @file rotating_buffer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include <cstddef>

namespace numint::detail
{

template <class T, std::size_t N>
class rotating_buffer {
public:
    using value_type                 = T;
    const static std::size_t dimension = N;

    rotating_buffer()
        : m_first(0)
    {
    }

    std::size_t size() const
    {
        return dimension;
    }

    value_type &operator[](std::size_t i)
    {
        return m_data[get_index(i)];
    }

    const value_type &operator[](std::size_t i) const
    {
        return m_data[get_index(i)];
    }

    void rotate()
    {
        if (m_first == 0)
            m_first = dimension - 1;
        else
            --m_first;
    }

protected:
    value_type m_data[N];

private:
    std::size_t get_index(std::size_t i) const
    {
        return ((i + m_first) % dimension);
    }

    std::size_t m_first;
};

} // namespace numint::detail
