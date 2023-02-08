# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
cmake_minimum_required(VERSION 3.1...3.18)

# Set the project name.
project(stunning_chainsaw CXX)

# Set the default build type to Debug.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

# -----------------------------------------------------------------------------
# OPTIONS
# -----------------------------------------------------------------------------

option(STN_CHN_BUILD_EXAMPLES "Build examples" OFF)
option(STN_CHN_STRICT_WARNINGS "Enable strict compiler warnings" ON)
option(STN_CHN_WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)
option(STN_CHN_ENABLE_PLOT "Enable matplot++ plots for examples" OFF)

# -----------------------------------------------------------------------------
# Set the compilation flags.
# -----------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings that suggest using MSVC-specific safe functions
    set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} -D_CRT_SECURE_NO_WARNINGS)

    if(STN_CHN_WARNINGS_AS_ERRORS)
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /WX)
    endif(STN_CHN_WARNINGS_AS_ERRORS)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(STN_CHN_WARNINGS_AS_ERRORS)
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} -Werror)
    endif(STN_CHN_WARNINGS_AS_ERRORS)
endif()

if(STN_CHN_STRICT_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # Mark system headers as external for MSVC explicitly
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /experimental:external)
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /external:I ${CMAKE_BINARY_DIR})
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /external:anglebrackets)
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /external:W0)
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} /W4)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS} -Wall -Wextra -Wconversion -pedantic)
    endif()
endif(STN_CHN_STRICT_WARNINGS)

# -----------------------------------------------------------------------------
# LIBRARY
# -----------------------------------------------------------------------------

# Add the C++ library.
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
# Inlcude header directories.
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/include)
# Set compilation flags.
target_compile_options(${PROJECT_NAME} INTERFACE ${STN_CHN_COMPILE_OPTIONS})

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------

if(STN_CHN_BUILD_EXAMPLES)

    # = FETCH =================================================================

    # We need this in order to import external projects.
    include(FetchContent)

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
            FETCHCONTENT_QUIET FETCHCONTENT_BASE_DIR FETCHCONTENT_FULLY_DISCONNECTED FETCHCONTENT_UPDATES_DISCONNECTED
            FETCHCONTENT_UPDATES_DISCONNECTED_STOPWATCH FETCHCONTENT_SOURCE_DIR_STOPWATCH
        )
    endif(NOT stopwatch_POPULATED)

    # = TARGETS ===============================================================

    # Add the lotka example.
    add_executable(lotka example/lotka.cpp)
    # Inlcude header directories.
    target_include_directories(lotka PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(lotka PUBLIC ${STN_CHN_COMPILE_OPTIONS})
    # Set linking flags.
    target_link_libraries(lotka PUBLIC ${PROJECT_NAME})

    # Add the dcmotor example.
    add_executable(dcmotor example/dcmotor.cpp)
    # Inlcude header directories.
    target_include_directories(dcmotor PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(dcmotor PUBLIC ${STN_CHN_COMPILE_OPTIONS})
    # Set linking flags.
    target_link_libraries(dcmotor PUBLIC ${PROJECT_NAME})

    # Add the dcmotor_v2 example.
    add_executable(dcmotor_v2 example/dcmotor_v2.cpp)
    # Inlcude header directories.
    target_include_directories(dcmotor_v2 PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(dcmotor_v2 PUBLIC ${STN_CHN_COMPILE_OPTIONS})
    # Set linking flags.
    target_link_libraries(dcmotor_v2 PUBLIC ${PROJECT_NAME})

    # Add the dcmotor example.
    add_executable(tandem_dc_motors example/tandem_dc_motors.cpp)
    # Inlcude header directories.
    target_include_directories(tandem_dc_motors PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(tandem_dc_motors PUBLIC ${STN_CHN_COMPILE_OPTIONS})
    # Set linking flags.
    target_link_libraries(tandem_dc_motors PUBLIC ${PROJECT_NAME})

    # Add the dcmotor example.
    add_executable(bounching_ball example/bounching_ball.cpp)
    # Inlcude header directories.
    target_include_directories(bounching_ball PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(bounching_ball PUBLIC ${STN_CHN_COMPILE_OPTIONS})
    # Set linking flags.
    target_link_libraries(bounching_ball PUBLIC ${PROJECT_NAME})

    # Add matplot++ if required.
    if(STN_CHN_ENABLE_PLOT)
        # Find matplot++ library.
        find_package(Matplot++ CONFIG REQUIRED)

        # Link libraries.
        target_link_libraries(lotka PUBLIC Matplot++::matplot)
        target_link_libraries(dcmotor PUBLIC Matplot++::matplot)
        target_link_libraries(dcmotor_v2 PUBLIC Matplot++::matplot)
        target_link_libraries(tandem_dc_motors PUBLIC Matplot++::matplot)
        target_link_libraries(bounching_ball PUBLIC Matplot++::matplot)

        # Add definitions.
        target_compile_definitions(lotka PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(dcmotor PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(dcmotor_v2 PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(tandem_dc_motors PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(bounching_ball PUBLIC SC_ENABLE_PLOT)
    endif(STN_CHN_ENABLE_PLOT)
endif(STN_CHN_BUILD_EXAMPLES)

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------
find_package(Doxygen)

if(DOXYGEN_FOUND)

    # = FETCH =================================================================

    # Include fetch content.
    include(FetchContent)
    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(
        doxygenawesome
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
        GIT_TAG 4cd62308d825fe0396d2f66ffbab45d0e247724c # 2.0.3
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
            FETCHCONTENT_QUIET FETCHCONTENT_BASE_DIR FETCHCONTENT_FULLY_DISCONNECTED FETCHCONTENT_UPDATES_DISCONNECTED
            FETCHCONTENT_UPDATES_DISCONNECTED_DOXYGENAWESOME FETCHCONTENT_SOURCE_DIR_DOXYGENAWESOME
        )
    endif()

    # = TARGETS ===============================================================

    set(DOXYGEN_PROJECT_NAME "Integration Library")
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
        ${PROJECT_SOURCE_DIR}/include/solver/solver.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/detail/it_algebra.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/detail/less_with_sign.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/detail/observer.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/detail/rotating_buffer.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/detail/type_traits.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/stepper/stepper_adaptive.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/stepper/stepper_euler.hpp
        ${PROJECT_SOURCE_DIR}/include/solver/stepper/stepper_rk4.hpp
    )
endif()
