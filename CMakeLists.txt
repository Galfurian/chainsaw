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
option(ENABLE_PLOT "Enable matplot++ plots for examples" OFF)
option(BUILD_EXAMPLES "Build examples" ON)

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------

# We want doxygen for the documentation.
find_package(Doxygen)

# Find clang-tidy if available.
find_program(CLANG_TIDY_EXE NAMES clang-tidy)

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
    # Disable warnings for MSVC-specific "safe" functions like strcpy_s, etc.,
    # which are not portable and may clutter warning logs.
    target_compile_definitions(${PROJECT_NAME} INTERFACE _CRT_SECURE_NO_WARNINGS)

    # Disable warning C4702: unreachable code.
    add_compile_options(/wd4702)

    if(WARNINGS_AS_ERRORS)
        # Treat all warnings as errors to enforce stricter code quality.
        target_compile_options(${PROJECT_NAME} INTERFACE /WX)
    endif()

    if(STRICT_WARNINGS)
        # Enable external header management to suppress warnings in system and
        # external headers, making it easier to focus on project-specific issues.
        target_compile_options(${PROJECT_NAME} INTERFACE /experimental:external)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:I ${CMAKE_BINARY_DIR})
        target_compile_options(${PROJECT_NAME} INTERFACE /external:anglebrackets)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:W0)

        # Use a high warning level to catch as many potential issues as possible.
        target_compile_options(${PROJECT_NAME} INTERFACE /W4)

        # Enforce standards-compliant behavior to avoid relying on MSVC-specific extensions.
        target_compile_options(${PROJECT_NAME} INTERFACE /permissive-)
    endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(WARNINGS_AS_ERRORS)
        # Treat all warnings as errors to enforce stricter code quality.
        target_compile_options(${PROJECT_NAME} INTERFACE -Werror)
    endif()

    if(STRICT_WARNINGS)
        # Enable a broad set of warnings to catch common and subtle issues:
        target_compile_options(${PROJECT_NAME} INTERFACE
            -Wall                # Enable most general-purpose warnings.
            -Wextra              # Enable extra warnings not included in -Wall.
            -Wconversion         # Warn about implicit type conversions that may lose data.
            -pedantic            # Enforce strict compliance with the C++ standard.
            -Wshadow             # Warn about variable shadowing, which can cause subtle bugs.
            -Wnon-virtual-dtor   # Warn when a class with virtual functions lacks a virtual destructor.
            -Wnull-dereference   # Warn about potential null pointer dereferences.
            -Wformat=2           # Enable strict checks for printf/scanf format strings.
            -Woverloaded-virtual # Warn when a derived class function hides a base class virtual function.
            -Wfloat-equal        # Warn about direct comparisons of floating-point values, which can be imprecise.
        )

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
# CODE ANALYSIS
# -----------------------------------------------------------------------------

if(CLANG_TIDY_EXE)
    file(GLOB_RECURSE ALL_PROJECT_FILES
        "${PROJECT_SOURCE_DIR}/include/**/*.hpp"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/**/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
    )
    add_custom_target(
        ${PROJECT_NAME}_clang_tidy
        COMMAND ${CLANG_TIDY_EXE}
        --p=${CMAKE_BINARY_DIR}
        ${ALL_PROJECT_FILES}
        COMMENT "Running clang-tidy"
        VERBATIM
    )
    add_custom_target(
        ${PROJECT_NAME}_clang_tidy_fix
        COMMAND ${CLANG_TIDY_EXE}
        --fix
        --p=${CMAKE_BINARY_DIR}
        ${ALL_PROJECT_FILES}
        COMMENT "Running clang-tidy-fix"
        VERBATIM
    )
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
    set(DOXYGEN_WARN_FORMAT "$file:$line:1: $text")
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
    set(DOXYGEN_EXCLUDE_PATTERNS "${PROJECT_SOURCE_DIR}/tests/*" "${PROJECT_SOURCE_DIR}/examples/*")

    # Add Doxygen documentation target.
    file(GLOB_RECURSE ALL_PROJECT_FILES
        "${PROJECT_SOURCE_DIR}/include/**/*.hpp"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/**/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
    )
    doxygen_add_docs(
        numint_documentation
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/LICENSE.md
        ${ALL_PROJECT_FILES}
        COMMENT "Generating Doxygen documentation for the library"
    )

endif()
