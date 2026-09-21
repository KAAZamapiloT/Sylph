cmake_minimum_required(VERSION 3.24)

project(
    Syplh
    VERSION 0.1.0
    DESCRIPTION "Spherical Harmonics and Spherical Grid Research Library"
    LANGUAGES C CXX
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
    ${CMAKE_SOURCE_DIR}/external/eigen
    ${CMAKE_BINARY_DIR}/eigen
    EXCLUDE_FROM_ALL
)

# ============================================================
# SDL3
# Prebuilt Windows x64 package
# ============================================================

set(SDL3_ROOT
    ${CMAKE_SOURCE_DIR}/external/SDL3-3.4.0
)

set(SDL3_INCLUDE_DIR
    ${SDL3_ROOT}/include
)

set(SDL3_LIBRARY
    ${SDL3_ROOT}/lib/x64/SDL3.lib
)

set(SDL3_DLL
    ${SDL3_ROOT}/lib/x64/SDL3.dll
)

add_library(SDL3::SDL3 SHARED IMPORTED GLOBAL)

set_target_properties(
    SDL3::SDL3
    PROPERTIES

        INTERFACE_INCLUDE_DIRECTORIES
            "${SDL3_INCLUDE_DIR}"

        IMPORTED_IMPLIB
            "${SDL3_LIBRARY}"

        IMPORTED_LOCATION
            "${SDL3_DLL}"
)

# ============================================================
# OpenGL
# ============================================================

find_package(OpenGL REQUIRED)

# ============================================================
# GLAD
# ============================================================

add_library(glad STATIC
    ${CMAKE_SOURCE_DIR}/external/glad/src/gl.c
)

target_include_directories(
    glad
    PUBLIC
        ${CMAKE_SOURCE_DIR}/external/glad/include
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
# Syplh Renderer
# ============================================================

add_library(SyplhRenderer STATIC

    ${CMAKE_SOURCE_DIR}/src/render/camera.cpp
    ${CMAKE_SOURCE_DIR}/src/render/transform.cpp
    ${CMAKE_SOURCE_DIR}/src/render/mesh.cpp
    ${CMAKE_SOURCE_DIR}/src/render/shader.cpp
    ${CMAKE_SOURCE_DIR}/src/render/material.cpp
    ${CMAKE_SOURCE_DIR}/src/render/renderer.cpp
)

target_include_directories(
    SyplhRenderer
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(
    SyplhRenderer
    PUBLIC
        Syplh
        SDL3::SDL3
        glad
        OpenGL::GL
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
        SyplhRenderer
)

# ============================================================
# Warnings
# ============================================================

if (MSVC)

    target_compile_options(
        SyplhRenderer
        PRIVATE
            /W4
            /permissive-
    )

    target_compile_options(
        syplh_demo
        PRIVATE
            /W4
            /permissive-
    )

else()

    target_compile_options(
        SyplhRenderer
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )

    target_compile_options(
        syplh_demo
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )

endif()

# ============================================================
# Copy SDL3.dll beside executable
# ============================================================

add_custom_command(
    TARGET syplh_demo
    POST_BUILD

    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "${SDL3_DLL}"
        "$<TARGET_FILE_DIR:syplh_demo>"
)