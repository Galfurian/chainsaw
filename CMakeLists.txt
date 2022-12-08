# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
cmake_minimum_required(VERSION 3.1...3.18)

# Set the project name.
project(stunning-chainsaw CXX)

# Set C++ Standard.
set(CMAKE_CXX_STANDARD 17)

# -----------------------------------------------------------------------------
# Add options.
# -----------------------------------------------------------------------------
option(STN_CHN_BUILD_EXAMPLES "Build stunning-chainsaw examples" OFF)
option(STN_CHN_STRICT_WARNINGS "Enable strict compiler warnings for stunning-chainsaw" ON)
option(STN_CHN_WARNINGS_AS_ERRORS "Treat all warnings as errors for stunning-chainsaw" OFF)
option(STN_CHN_ENABLE_PLOT "Enable matplot++ plots for stunning-chainsaw examples" OFF)

# -----------------------------------------------------------------------------
# Set the compilation flags.
# -----------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings that suggest using MSVC-specific safe functions
    set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} -D_CRT_SECURE_NO_WARNINGS")

    if(STN_CHN_WARNINGS_AS_ERRORS)
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /WX")
    endif(STN_CHN_WARNINGS_AS_ERRORS)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(STN_CHN_WARNINGS_AS_ERRORS)
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} -Werror")
    endif(STN_CHN_WARNINGS_AS_ERRORS)
endif()

if(STN_CHN_STRICT_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # Mark system headers as external for MSVC explicitly
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /experimental:external")
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /external:I ${CMAKE_BINARY_DIR}")
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /external:anglebrackets")
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /external:W0")
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} /W4")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(STN_CHN_COMPILE_OPTIONS "${STN_CHN_COMPILE_OPTIONS} -Wall -Wextra -Wconversion -pedantic")
    endif()
endif(STN_CHN_STRICT_WARNINGS)

string(REPLACE " " ";" STN_CHN_COMPILE_OPTIONS ${STN_CHN_COMPILE_OPTIONS})

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------
find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME "Integration Library")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${PROJECT_SOURCE_DIR}/README.md)
    set(DOXYGEN_SHOW_INCLUDE_FILES NO)
    set(DOXYGEN_GENERATE_TREEVIEW NO)
    set(DOXYGEN_HTML_STYLESHEET ${PROJECT_SOURCE_DIR}/docs/doxygen.css)
    doxygen_add_docs(
        ${PROJECT_NAME}_documentation
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

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------
if(STN_CHN_BUILD_EXAMPLES)
    # We need this in order to import external projects.
    include(FetchContent)

    # Retrieve the stopwatch.
    message(STATUS "Retrieving `stopwatch`...")
    FetchContent_Declare(stopwatch
        GIT_REPOSITORY "https://github.com/Galfurian/stopwatch.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    FetchContent_GetProperties(stopwatch)

    if(NOT stopwatch_POPULATED)
        FetchContent_MakeAvailable(stopwatch)
    endif(NOT stopwatch_POPULATED)

    # Add the lotka example.
    add_executable(lotka example/lotka.cpp)
    # Inlcude header directories.
    target_include_directories(lotka PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(lotka PUBLIC ${STN_CHN_COMPILE_OPTIONS})

    # Add the dcmotor example.
    add_executable(dcmotor example/dcmotor.cpp)
    # Inlcude header directories.
    target_include_directories(dcmotor PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(dcmotor PUBLIC ${STN_CHN_COMPILE_OPTIONS})

    # Add the dcmotor_v2 example.
    add_executable(dcmotor_v2 example/dcmotor_v2.cpp)
    # Inlcude header directories.
    target_include_directories(dcmotor_v2 PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(dcmotor_v2 PUBLIC ${STN_CHN_COMPILE_OPTIONS})

    # Add the dcmotor example.
    add_executable(tandem_dc_motors example/tandem_dc_motors.cpp)
    # Inlcude header directories.
    target_include_directories(tandem_dc_motors PUBLIC ${PROJECT_SOURCE_DIR}/include ${stopwatch_SOURCE_DIR}/include)
    # Set compilation flags.
    target_compile_options(tandem_dc_motors PUBLIC ${STN_CHN_COMPILE_OPTIONS})

    # Add matplot++ if required.
    if(STN_CHN_ENABLE_PLOT)
        # Find matplot++ library.
        find_package(Matplot++ CONFIG REQUIRED)

        # Link libraries.
        target_link_libraries(lotka PUBLIC Matplot++::matplot)
        target_link_libraries(dcmotor PUBLIC Matplot++::matplot)
        target_link_libraries(dcmotor_v2 PUBLIC Matplot++::matplot)
        target_link_libraries(tandem_dc_motors PUBLIC Matplot++::matplot)

        # Add definitions.
        target_compile_definitions(lotka PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(dcmotor PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(dcmotor_v2 PUBLIC SC_ENABLE_PLOT)
        target_compile_definitions(tandem_dc_motors PUBLIC SC_ENABLE_PLOT)
    endif(STN_CHN_ENABLE_PLOT)
endif(STN_CHN_BUILD_EXAMPLES)
