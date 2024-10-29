# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
cmake_minimum_required(VERSION 3.1...3.18)

# Set the project name.
project(chainsaw CXX)

# Set the default build type to Debug.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

# -----------------------------------------------------------------------------
# ENABLE FETCH CONTENT
# -----------------------------------------------------------------------------

# We need this in order to import external projects.
include(FetchContent)
# Hide fetchcontent variables.
mark_as_advanced(FORCE
    FETCHCONTENT_QUIET
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_UPDATES_DISCONNECTED
)

# -----------------------------------------------------------------------------
# OPTIONS
# -----------------------------------------------------------------------------

option(STRICT_WARNINGS "Enable strict compiler warnings" ON)
option(WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(ENABLE_PLOT "Enable matplot++ plots for examples" OFF)

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------

# We want doxygen for the documentation.
find_package(Doxygen)

# -----------------------------------------------------------------------------
# LIBRARY
# -----------------------------------------------------------------------------

# Add the C++ library.
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
# Inlcude header directories.
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/include)
# Set the library to use c++-17
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

# -----------------------------------------------------------------------------
# Set the compilation flags.
# -----------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings that suggest using MSVC-specific safe functions
    target_compile_definitions(${PROJECT_NAME} INTERFACE _CRT_SECURE_NO_WARNINGS)
    if(WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} INTERFACE /WX)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} INTERFACE -Werror)
    endif()
endif()

if(STRICT_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # Mark system headers as external for MSVC explicitly
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory
        target_compile_options(${PROJECT_NAME} INTERFACE /experimental:external)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:I ${CMAKE_BINARY_DIR})
        target_compile_options(${PROJECT_NAME} INTERFACE /external:anglebrackets)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:W0)

        target_compile_options(${PROJECT_NAME} INTERFACE /W4)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${PROJECT_NAME} INTERFACE -Wall -Wextra -Wconversion -pedantic)
    endif()
endif()

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------

if(BUILD_EXAMPLES)
        
    # = FETCH =================================================================
    # Retrieve the stopwatch.
    FetchContent_Declare(stopwatch
        GIT_REPOSITORY "https://github.com/Galfurian/stopwatch.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    FetchContent_GetProperties(stopwatch)
    if(NOT stopwatch_POPULATED)
        message(STATUS "Retrieving `stopwatch`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(stopwatch)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_STOPWATCH
            FETCHCONTENT_SOURCE_DIR_STOPWATCH
        )
    endif(NOT stopwatch_POPULATED)

    # = TARGETS ===============================================================

    # Add the example.
    add_executable(${PROJECT_NAME}_lotka ${PROJECT_SOURCE_DIR}/examples/lotka.cpp)
    target_include_directories(${PROJECT_NAME}_lotka PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_lotka PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_dcmotor ${PROJECT_SOURCE_DIR}/examples/dcmotor.cpp)
    target_include_directories(${PROJECT_NAME}_dcmotor PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_dcmotor PUBLIC ${PROJECT_NAME})

    # Add the dcmotor_v2 example.
    add_executable(${PROJECT_NAME}_dcmotor_v2 ${PROJECT_SOURCE_DIR}/examples/dcmotor_v2.cpp)
    target_include_directories(${PROJECT_NAME}_dcmotor_v2 PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_dcmotor_v2 PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_tandem_dc_motors ${PROJECT_SOURCE_DIR}/examples/tandem_dc_motors.cpp)
    target_include_directories(${PROJECT_NAME}_tandem_dc_motors PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_tandem_dc_motors PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_bounching_ball ${PROJECT_SOURCE_DIR}/examples/bounching_ball.cpp)
    target_include_directories(${PROJECT_NAME}_bounching_ball PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_bounching_ball PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_spring_mass_damper ${PROJECT_SOURCE_DIR}/examples/spring_mass_damper.cpp)
    target_include_directories(${PROJECT_NAME}_spring_mass_damper PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_spring_mass_damper PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_pendulum ${PROJECT_SOURCE_DIR}/examples/pendulum.cpp)
    target_include_directories(${PROJECT_NAME}_pendulum PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_accellerometer ${PROJECT_SOURCE_DIR}/examples/accellerometer.cpp)
    target_include_directories(${PROJECT_NAME}_accellerometer PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_accellerometer PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_compare_fixed ${PROJECT_SOURCE_DIR}/examples/compare_fixed.cpp)
    target_include_directories(${PROJECT_NAME}_compare_fixed PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_compare_fixed PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_compare_adaptive ${PROJECT_SOURCE_DIR}/examples/compare_adaptive.cpp)
    target_include_directories(${PROJECT_NAME}_compare_adaptive PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_compare_adaptive PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_robot_arm ${PROJECT_SOURCE_DIR}/examples/robot_arm.cpp)
    target_include_directories(${PROJECT_NAME}_robot_arm PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_robot_arm PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_multi_mode ${PROJECT_SOURCE_DIR}/examples/multi_mode.cpp)
    target_include_directories(${PROJECT_NAME}_multi_mode PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_multi_mode PUBLIC ${PROJECT_NAME})

    # Add matplot++ if required.
    if(ENABLE_PLOT)
        # Find matplot++ library.
        find_package(Matplot++ CONFIG REQUIRED)

        # Link libraries.
        target_link_libraries(${PROJECT_NAME}_lotka PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_dcmotor PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_dcmotor_v2 PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_tandem_dc_motors PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_bounching_ball PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_spring_mass_damper PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_accellerometer PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_compare_fixed PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_compare_adaptive PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_robot_arm PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_multi_mode PUBLIC Matplot++::matplot)

        # Add definitions.
        target_compile_definitions(${PROJECT_NAME}_lotka PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_dcmotor PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_dcmotor_v2 PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_tandem_dc_motors PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_bounching_ball PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_spring_mass_damper PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_pendulum PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_accellerometer PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_compare_fixed PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_compare_adaptive PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_robot_arm PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_multi_mode PUBLIC SC_ENABLE_PLOT)
    endif()
endif()

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------

if(DOXYGEN_FOUND)

    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(
        doxygenawesome
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
        GIT_TAG main
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(doxygenawesome)
    # If not populated, make the content available.
    if(NOT doxygenawesome_POPULATED)
        message(STATUS "Retrieving `doxygen-awesome-css`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(doxygenawesome)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_DOXYGENAWESOME
            FETCHCONTENT_SOURCE_DIR_DOXYGENAWESOME
        )
    endif()

    # = CUSTOMIZATION =========================================================
    set(DOXYGEN_WARN_FORMAT "$file:$line: $text")
    set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${PROJECT_SOURCE_DIR}/README.md)
    set(DOXYGEN_SHOW_INCLUDE_FILES NO)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_HTML_HEADER ${doxygenawesome_SOURCE_DIR}/doxygen-custom/header.html)
    set(DOXYGEN_HTML_EXTRA_STYLESHEET ${doxygenawesome_SOURCE_DIR}/doxygen-awesome.css)
    set(DOXYGEN_HTML_EXTRA_FILES
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
    )
    doxygen_add_docs(
        ${PROJECT_NAME}_documentation
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/include/chainsaw/solver.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/detail/it_algebra.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/detail/less_with_sign.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/detail/observer.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/detail/rotating_buffer.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/detail/type_traits.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/stepper/stepper_adaptive.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/stepper/stepper_euler.hpp
        ${PROJECT_SOURCE_DIR}/include/chainsaw/stepper/stepper_rk4.hpp
    )
endif()
