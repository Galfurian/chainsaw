/// @file rotating_buffer.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Simplification of the code available at:
///     https://github.com/headmyshoulder/odeint-v2

#pragma once

#include <array>
#include <cstddef>

namespace numint::detail
{

/// @brief A rotating buffer.
/// @tparam T the type of the elements.
/// @tparam N the number of elements.
template <class T, std::size_t N>
class rotating_buffer
{
public:
    /// @brief The type of the elements.
    using value_type                   = T;
    /// @brief The number of elements.
    const static std::size_t dimension = N;

    /// @brief Construct a new rotating buffer.
    rotating_buffer() = default;

    /// @brief Get the size of the buffer.
    /// @return The size of the buffer.
    auto size() const -> std::size_t { return dimension; }

    /// @brief Get the element at the given index.
    /// @param index the index.
    /// @return The element at the given index.
    auto operator[](std::size_t index) -> value_type & { return m_data[get_index(index)]; }

    /// @brief Get the element at the given index.
    /// @param index the index.
    /// @return The element at the given index.
    auto operator[](std::size_t index) const -> const value_type & { return m_data[get_index(index)]; }

    /// @brief Rotate the buffer.
    void rotate()
    {
        if (m_first == 0) {
            m_first = dimension - 1;
        } else {
            --m_first;
        }
    }

private:
    /// @brief The data.
    std::array<value_type, N> m_data;

    /// @brief Returns the correct index of the element at the given index.
    /// @param index the index.
    /// @return The correct index of the element at the given index.
    auto get_index(std::size_t index) const -> std::size_t { return ((index + m_first) % dimension); }

    /// @brief The first element.
    std::size_t m_first{0};
};

} // namespace numint::detail
