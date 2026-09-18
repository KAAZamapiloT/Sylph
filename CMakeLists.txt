cmake_minimum_required(VERSION 3.24)

project(
    Syplh
    VERSION 0.1.0
    DESCRIPTION "Spherical Harmonics and Spherical Grid Research Library"
    LANGUAGES CXX
)

# ============================================================
# C++23
# ============================================================

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================
# Eigen
# ============================================================

add_subdirectory(
    ${CMAKE_SOURCE_DIR}/External/eigen
    ${CMAKE_BINARY_DIR}/eigen
    EXCLUDE_FROM_ALL
)

# ============================================================
# Syplh - Header-only library
# ============================================================

add_library(Syplh INTERFACE)

target_compile_features(
    Syplh
    INTERFACE
        cxx_std_23
)

target_include_directories(
    Syplh
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(
    Syplh
    INTERFACE
        Eigen3::Eigen
)

# ============================================================
# Development executable
# ============================================================

add_executable(
    syplh_demo
    examples/main.cpp
)

target_link_libraries(
    syplh_demo
    PRIVATE
        Syplh
)

# ============================================================
# Warnings
# ============================================================

if (MSVC)
    target_compile_options(
        syplh_demo
        PRIVATE
            /W4
            /permissive-
    )
else()
    target_compile_options(
        syplh_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()
