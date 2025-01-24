/// @file rotating_buffer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include <cstddef>

namespace numint::detail
{

/// @brief A rotating buffer.
/// @tparam T the type of the elements.
/// @tparam N the number of elements.
template <class T, std::size_t N>
class rotating_buffer {
public:
    /// @brief The type of the elements.
    using value_type                   = T;
    /// @brief The number of elements.
    const static std::size_t dimension = N;

    /// @brief Construct a new rotating buffer.
    rotating_buffer()
        : m_first(0)
    {
    }

    /// @brief Get the size of the buffer.
    /// @return The size of the buffer.
    std::size_t size() const
    {
        return dimension;
    }

    /// @brief Get the element at the given index.
    /// @param i the index.
    /// @return The element at the given index.
    value_type &operator[](std::size_t i)
    {
        return m_data[get_index(i)];
    }

    /// @brief Get the element at the given index.
    /// @param i the index.
    /// @return The element at the given index.
    const value_type &operator[](std::size_t i) const
    {
        return m_data[get_index(i)];
    }

    /// @brief Rotate the buffer.
    void rotate()
    {
        if (m_first == 0)
            m_first = dimension - 1;
        else
            --m_first;
    }

protected:
    /// @brief The data.
    value_type m_data[N];

private:
    /// @brief Returns the correct index of the element at the given index.
    /// @param i the index.
    /// @return The correct index of the element at the given index.
    std::size_t get_index(std::size_t i) const
    {
        return ((i + m_first) % dimension);
    }

    /// @brief The first element.
    std::size_t m_first;
};

} // namespace numint::detail
