# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
cmake_minimum_required(VERSION 3.1...3.18)

# Set the project name.
project(numint CXX)

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
    # Retrieve the timelib.
    FetchContent_Declare(timelib
        GIT_REPOSITORY "https://github.com/Galfurian/timelib.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    FetchContent_GetProperties(timelib)
    if(NOT timelib_POPULATED)
        message(STATUS "Retrieving `timelib`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(timelib)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_TIMELIB
            FETCHCONTENT_SOURCE_DIR_TIMELIB
        )
    endif(NOT timelib_POPULATED)

    # = TARGETS ===============================================================

    # Add the example.
    add_executable(${PROJECT_NAME}_lotka ${PROJECT_SOURCE_DIR}/examples/lotka.cpp)
    target_include_directories(${PROJECT_NAME}_lotka PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_lotka PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_dcmotor ${PROJECT_SOURCE_DIR}/examples/dcmotor.cpp)
    target_include_directories(${PROJECT_NAME}_dcmotor PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_dcmotor PUBLIC ${PROJECT_NAME})

    # Add the dcmotor_v2 example.
    add_executable(${PROJECT_NAME}_dcmotor_v2 ${PROJECT_SOURCE_DIR}/examples/dcmotor_v2.cpp)
    target_include_directories(${PROJECT_NAME}_dcmotor_v2 PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_dcmotor_v2 PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_tandem_dc_motors ${PROJECT_SOURCE_DIR}/examples/tandem_dc_motors.cpp)
    target_include_directories(${PROJECT_NAME}_tandem_dc_motors PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_tandem_dc_motors PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_bounching_ball ${PROJECT_SOURCE_DIR}/examples/bounching_ball.cpp)
    target_include_directories(${PROJECT_NAME}_bounching_ball PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_bounching_ball PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_spring_mass_damper ${PROJECT_SOURCE_DIR}/examples/spring_mass_damper.cpp)
    target_include_directories(${PROJECT_NAME}_spring_mass_damper PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_spring_mass_damper PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_pendulum ${PROJECT_SOURCE_DIR}/examples/pendulum.cpp)
    target_include_directories(${PROJECT_NAME}_pendulum PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC ${PROJECT_NAME})

    # Add the example.
    add_executable(${PROJECT_NAME}_compare_fixed ${PROJECT_SOURCE_DIR}/examples/compare_fixed.cpp)
    target_include_directories(${PROJECT_NAME}_compare_fixed PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_compare_fixed PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_compare_adaptive ${PROJECT_SOURCE_DIR}/examples/compare_adaptive.cpp)
    target_include_directories(${PROJECT_NAME}_compare_adaptive PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_compare_adaptive PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_robot_arm ${PROJECT_SOURCE_DIR}/examples/robot_arm.cpp)
    target_include_directories(${PROJECT_NAME}_robot_arm PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_robot_arm PUBLIC ${PROJECT_NAME})
    
    # Add the example.
    add_executable(${PROJECT_NAME}_multi_mode ${PROJECT_SOURCE_DIR}/examples/multi_mode.cpp)
    target_include_directories(${PROJECT_NAME}_multi_mode PUBLIC ${PROJECT_SOURCE_DIR}/include ${timelib_SOURCE_DIR}/include)
    target_link_libraries(${PROJECT_NAME}_multi_mode PUBLIC ${PROJECT_NAME})

    # Add matplot++ if required.
    if(ENABLE_PLOT)
    
        FetchContent_Declare(
            gpcpp
            GIT_REPOSITORY "https://github.com/Galfurian/gpcpp.git"
            GIT_SHALLOW    TRUE
            GIT_PROGRESS   TRUE
            GIT_TAG        main
        )
        FetchContent_GetProperties(gpcpp)
        if(NOT gpcpp_POPULATED)
            message(STATUS "Retrieving `gpcpp`...")
            FetchContent_MakeAvailable(gpcpp)
            # Hide fetchcontent variables.
            mark_as_advanced(FORCE FETCHCONTENT_UPDATES_DISCONNECTED_GPCPP FETCHCONTENT_SOURCE_DIR_GPCPP)
        endif()


        # Link libraries.
        target_link_libraries(${PROJECT_NAME}_lotka PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_dcmotor PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_dcmotor_v2 PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_tandem_dc_motors PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_bounching_ball PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_spring_mass_damper PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_compare_fixed PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_compare_adaptive PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_robot_arm PUBLIC gpcpp)
        target_link_libraries(${PROJECT_NAME}_multi_mode PUBLIC gpcpp)

        # Add definitions.
        target_compile_definitions(${PROJECT_NAME}_lotka PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_dcmotor PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_dcmotor_v2 PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_tandem_dc_motors PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_bounching_ball PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_spring_mass_damper PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_pendulum PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_compare_fixed PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_compare_adaptive PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_robot_arm PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_multi_mode PUBLIC ENABLE_PLOT)
    endif()
endif()

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------

if(DOXYGEN_FOUND)
    
    # Read the file with the version.
    file(READ ${PROJECT_SOURCE_DIR}/include/numint/solver.hpp version_file)
    # Extract the version.
    string(REGEX MATCH "NUMINT_MAJOR_VERSION ([0-9]*)" _ ${version_file})
    set(NUMINT_MAJOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "NUMINT_MINOR_VERSION ([0-9]*)" _ ${version_file})
    set(NUMINT_MINOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "NUMINT_MICRO_VERSION ([0-9]*)" _ ${version_file})
    set(NUMINT_MICRO_VERSION ${CMAKE_MATCH_1})

    # Customization: Doxygen Configuration
    set(DOXYGEN_WARN_FORMAT "$file:$line: $text")
    set(DOXYGEN_PROJECT_NAME numint)
    set(DOXYGEN_PROJECT_BRIEF "A library for numerical integration")
    set(DOXYGEN_PROJECT_NUMBER "${NUMINT_MAJOR_VERSION}.${NUMINT_MINOR_VERSION}.${NUMINT_MICRO_VERSION}")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${PROJECT_SOURCE_DIR}/README.md)
    set(DOXYGEN_SHOW_INCLUDE_FILES NO)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_GENERATE_LATEX NO)
    set(DOXYGEN_GENERATE_MAN NO)  

    # FetchContent: Doxygen Awesome CSS
    FetchContent_Declare(doxygenawesome
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
        GIT_TAG main
    )
    FetchContent_MakeAvailable(doxygenawesome)

    # Hide FetchContent variables to avoid clutter in ccmake.
    mark_as_advanced(FORCE
        FETCHCONTENT_UPDATES_DISCONNECTED_DOXYGENAWESOME
        FETCHCONTENT_SOURCE_DIR_DOXYGENAWESOME
    )

    # Styling and UX enhancements using Doxygen Awesome
    set(DOXYGEN_HTML_HEADER ${doxygenawesome_SOURCE_DIR}/doxygen-custom/header.html)
    set(DOXYGEN_HTML_EXTRA_STYLESHEET ${doxygenawesome_SOURCE_DIR}/doxygen-awesome.css)
    set(DOXYGEN_HTML_EXTRA_FILES
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
    )

    # Set stricter warnings for better documentation quality
    set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
    set(DOXYGEN_WARN_IF_DOC_ERROR YES)
    set(DOXYGEN_WARN_NO_PARAMDOC YES)
    set(DOXYGEN_WARN_AS_ERROR YES) # Treat warnings as errors for CI

    # Exclude certain files or directories from documentation (if needed)
    set(DOXYGEN_EXCLUDE_PATTERNS
        "${PROJECT_SOURCE_DIR}/tests/*"
        "${PROJECT_SOURCE_DIR}/examples/*"
    )

    # Add Doxygen documentation target
    doxygen_add_docs(
        numint_documentation
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/include/numint/solver.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/detail/it_algebra.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/detail/less_with_sign.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/detail/observer.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/detail/rotating_buffer.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/detail/type_traits.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/stepper/stepper_adaptive.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/stepper/stepper_euler.hpp
        ${PROJECT_SOURCE_DIR}/include/numint/stepper/stepper_rk4.hpp
        COMMENT "Generating Doxygen documentation for the library"
    )

endif()
