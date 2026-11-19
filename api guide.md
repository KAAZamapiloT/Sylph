# Sylph: C++ API Reference & Developer Guide

A practical, hands-on developer guide and API reference for the **Sylph** Spherical Harmonics library.

---

## Table of Contents

1. [Quickstart & CMake Integration](#1-quickstart--cmake-integration)
2. [Core Class: `sylph::SHCoefficients`](#2-core-class-sylphshcoefficients)
3. [Core Class: `sylph::ZHCoefficients` (Zonal Harmonics)](#3-core-class-sylphzhcoefficients-zonal-harmonics)
4. [Transform Engine: `sylph::SHGridTransform`](#4-transform-engine-sylphshgridtransform)
5. [Spatial Grid: `sylph::Spherical_grid`](#5-spatial-grid-sylphspherical_grid)
6. [Low-Level Acceleration Caches](#6-low-level-acceleration-caches)
7. [Complete End-to-End Code Examples](#7-complete-end-to-end-code-examples)

---

## 1. Quickstart & CMake Integration

### Include Paths & Namespaces
Sylph is a modern C++20 header-first library. All core mathematical classes reside in the `sylph` namespace, while rendering extensions reside in the `render` namespace.

```cpp
#include "sylph/sh_coefficients.hpp"
#include "sylph/zonal_harmonics.hpp"
#include "sylph/sh_grid_transform.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/spherical_grid.hpp"

using namespace sylph;
```

### CMake Integration
Add Sylph to your `CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.20)
project(MySHApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add Sylph include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/Sylph/include
    ${CMAKE_CURRENT_SOURCE_DIR}/Sylph/include/types
)

# Optional: Enable AVX2 / OpenMP optimizations
if(MSVC)
    add_compile_options(/O2 /arch:AVX2 /fp:fast)
else()
    add_compile_options(-O3 -mavx2 -mfma -ffast-math)
endif()

add_executable(my_app main.cpp)
```

---

## 2. Core Class: `sylph::SHCoefficients`

The `SHCoefficients` class stores real Spherical Harmonic coefficients up to a given band-limit (order $N$).
- **Number of bands:** $l \in [0, N-1]$
- **Total stored coefficients:** $N^2$
- **Internal storage layout:** Flat 1D array indexed by:
  $$\text{index}(l, m) = l(l + 1) + m$$

### Key Member Functions

| Method | Return Type | Description |
|---|---|---|
| `SHCoefficients(int order)` | Constructor | Allocates an order-$N$ vector initialized to $0.0$. |
| `order()` | `int` | Returns the order $N$ (maximum degree is $N-1$). |
| `num_coefficients()` | `int` | Returns the total coefficient count ($N^2$). |
| `operator()(int l, int m)` | `double&` | Access coefficient for band $l$ and order $m \in [-l, l]$. |
| `operator()(int l, int m) const` | `double` | Read coefficient with boundary assertion. |
| `data()` | `double*` | Returns raw pointer to the contiguous double buffer. |
| `data() const` | `const double*` | Returns const raw pointer to the buffer. |
| `is_zonal(double eps = 1e-12)` | `bool` | Checks in $O(1)$ early-exit whether all $m \ne 0$ coefficients are zero. |

### Usage Example
```cpp
#include "sylph/sh_coefficients.hpp"
#include <iostream>

int main() {
    // Create an Order 4 SH representation (4^2 = 16 coefficients)
    sylph::SHCoefficients c(4);

    // Set band 0 (l=0, m=0): Monopole / Ambient light
    c(0, 0) = 1.0;

    // Set band 1 (l=1, m=1): Directional gradient along X-axis
    c(1, 1) = 0.5;

    // Set band 1 (l=1, m=-1): Directional gradient along Y-axis
    c(1, -1) = -0.25;

    std::cout << "Order: " << c.order() << "\n";
    std::cout << "Coefficients count: " << c.num_coefficients() << "\n";
    std::cout << "Is zonal: " << (c.is_zonal() ? "yes" : "no") << "\n";

    return 0;
}
```

---

## 3. Core Class: `sylph::ZHCoefficients` (Zonal Harmonics)

Zonal Harmonics represent functions that are rotationally symmetric around an axis. Since they have non-zero values **only for $m = 0$**, they store only $N$ numbers instead of $N^2$, enabling extreme $O(N^2)$ speedups.

### Key Member Functions

| Method | Return Type | Description |
|---|---|---|
| `ZHCoefficients(int order)` | Constructor | Allocates an order-$N$ Zonal Harmonic vector (size $N$). |
| `ZHCoefficients(span<const double>)` | Constructor | Constructs from an existing slice of $m=0$ coefficients. |
| `order()` | `int` | Returns the order $N$. |
| `operator[](int l)` | `double&` | Access degree $l \in [0, N-1]$ coefficient. |
| `data()` | `const double*` | Raw pointer to coefficients. |
| `rotate_to_sh(double theta, double phi)` | `SHCoefficients` | Rotates the axial ZH lobe to an arbitrary 3D direction in exact $O(N^2)$ time via the Addition Theorem. |
| `convolve(const SHCoefficients& F)` | `SHCoefficients` | Spherical convolution $(Z * F)$ in exact $O(N^2)$ time via the Funk–Hecke theorem. |

### Usage Example
```cpp
#include "sylph/zonal_harmonics.hpp"
#include "sylph/sh_coefficients.hpp"
#include <cmath>

int main() {
    const int order = 16;
    sylph::ZHCoefficients zh(order);

    // Define a clamped cosine lobe: z_l = ...
    zh[0] = 0.8862; // sqrt(pi)/2
    zh[1] = 1.0233; // sqrt(pi/3)
    zh[2] = 0.4954; // sqrt(5*pi)/8

    // 1. Rotate the lobe to point toward theta = 45 deg, phi = 90 deg (O(N^2))
    sylph::SHCoefficients rotated_lobe = zh.rotate_to_sh(M_PI / 4.0, M_PI / 2.0);

    // 2. Convolve with an environment lighting map in SH space (O(N^2))
    sylph::SHCoefficients env_light(order);
    env_light(0, 0) = 2.0;
    sylph::SHCoefficients convolved = zh.convolve(env_light);

    return 0;
}
```

---

## 4. Transform Engine: `sylph::SHGridTransform`

The core workhorse implementing the paper's Generalized Spherical Harmonic Product and Spherical Grid conversions.

### 4.1 Multiplying Two Functions (`product`)
```cpp
// Exact product of two functions truncated to output_order
SHCoefficients product(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order
) const;
```
- **Automated Dispatch:** If either `a` or `b` is detected as a Zonal Harmonic (`is_zonal()`), it automatically dispatches to the microsecond fused polar fast path (**$3.5\ \mu\text{s}$**).
- **Exactness:** Automatically computes exact quadrature degree $t = n_a + n_b + p - 2$ and configures the optimal Gauss–Legendre $\times$ FFT grid.

### 4.2 Generalized $K$-Way Product (`product`)
```cpp
SHCoefficients product(
    std::span<const SHCoefficients* const> inputs,
    int output_order
) const;
```
- Multiplies $k$ functions simultaneously: $G = \prod_{i=1}^k F_i$.
- **Asymmetric orders:** Each input pointer can point to an `SHCoefficients` of a completely different order $n_i$.
- Evaluates the product in a single combined grid pass, preventing repeated truncation error accumulation.

### 4.3 Approximate Mode for Extreme Speed (`product` with parameter `t`)
```cpp
SHCoefficients product(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order,
    int t
) const;
```
- Allows manually setting the quadrature resolution parameter $t$.
- Setting $t < n_G + p - 1$ trades higher-frequency precision for maximum real-time performance.

### 4.4 Product Integrals (`product_integral`)
```cpp
// Evaluates the integral of the product over the entire sphere
double product_integral(
    const SHCoefficients& a,
    const SHCoefficients& b
) const;

double product_integral(
    std::span<const SHCoefficients* const> inputs
) const;
```
- Evaluates $\int_{S^2} \prod_k F_k(\omega)\,d\Omega = \sqrt{4\pi}\,g_{0,0}$.
- For two inputs ($K = 2$), computes the exact analytical dot product $\mathbf{a} \cdot \mathbf{b}$ in nanoseconds.

---

## 5. Spatial Grid: `sylph::Spherical_grid`

Represents discrete samples on the 2-sphere sampled at Gauss–Legendre latitude nodes $\theta_i$ and uniform longitude nodes $\phi_j$.

```cpp
#include "sylph/spherical_grid.hpp"

// Create a grid with 16 latitude rings and 32 longitude samples
sylph::Spherical_grid grid(16, 32);

// Access grid properties
int nTheta = grid.theta_size();
int nPhi   = grid.phi_size();

// Read latitude and longitude angles
double lat = grid.theta(0); // colatitude in radians
double lon = grid.phi(0);   // longitude in radians

// Access pixel value at (i, j)
grid(0, 0) = 1.5;

// Pointwise multiplication in-place
grid.multiply_inplace(another_grid);
```

---

## 6. Low-Level Acceleration Caches

For performance-critical code requiring microsecond execution, Sylph exposes its internal thread-local static lookup caches:

### 6.1 Gauss–Legendre Cache (`gauss_legendre.hpp`)
```cpp
#include "sylph/gauss_legendre.hpp"

// O(1) thread-local lookup; zero Newton iterations
const auto& rule = sylph::gauss_legendre::get(16);

// Access precomputed nodes, weights, theta, and sin(theta)
double x0     = rule.nodes[0];
double w0     = rule.weights[0];
double theta0 = rule.theta[0];
double sin0   = rule.sin_theta[0];
```

### 6.2 Polar Basis Cache (`polar_basis_cache.hpp`)
```cpp
#include "sylph/polar_basis_cache.hpp"

// Precomputed Associated Legendre basis matrix in contiguous L1 cache layout
const auto& basis_cache = sylph::PolarBasisCache::get(nTheta, order);

// Get direct pointer to column values P~_l^m(x_i)
const double* p_col = basis_cache.column(m, k);
```

---

## 7. Complete End-to-End Code Examples

### Example 1: Multiplying Two General Order-16 Spherical Harmonics
```cpp
#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include <iostream>
#include <chrono>

int main() {
    const int order = 16;
    sylph::SHCoefficients lighting(order);
    sylph::SHCoefficients visibility(order);

    // Initialize sample coefficients
    lighting(0, 0) = 2.0;
    lighting(1, 0) = 0.5;
    visibility(0, 0) = 0.8;
    visibility(1, 0) = -0.3;

    sylph::SHGridTransform transform;

    auto start = std::chrono::high_resolution_clock::now();

    // Compute exact Order-16 product
    sylph::SHCoefficients result = transform.product(lighting, visibility, order);

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_us = std::chrono::duration<double, std::micro>(end - start).count();

    std::cout << "Product computed in: " << elapsed_us << " microseconds\n";
    std::cout << "Result monopole (0,0): " << result(0, 0) << "\n";

    return 0;
}
```

---

### Example 2: Ultra-Fast Zonal Harmonic × General SH Product
```cpp
#include "sylph/zonal_harmonics.hpp"
#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include <iostream>

int main() {
    const int order = 16;

    // 1. Create a directional sun lobe as a Zonal Harmonic
    sylph::ZHCoefficients sun_lobe(order);
    sun_lobe[0] = 3.0;
    sun_lobe[1] = 2.5;
    sun_lobe[2] = 1.0;

    // 2. General visibility field from an environment
    sylph::SHCoefficients env_visibility(order);
    env_visibility(0, 0) = 1.0;
    env_visibility(2, 2) = 0.4;

    sylph::SHGridTransform transform;

    // Automated dispatch detects sun_lobe as Zonal:
    // Runs in ~3.5 microseconds (fused polar quadrature, 0 FFTs, 0 2D grid)
    sylph::SHCoefficients lit_vis = transform.product(sun_lobe, env_visibility, order);

    std::cout << "ZH x SH Product successful. Band 0: " << lit_vis(0, 0) << "\n";
    return 0;
}
```

---

### Example 3: Generalized Triple Product with Mixed Input Orders
```cpp
#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include <array>
#include <iostream>

int main() {
    // Three inputs with completely different SH orders:
    // 1. Diffuse visibility: Order 4 (16 coeffs)
    // 2. Incident lighting: Order 16 (256 coeffs)
    // 3. Specular BRDF:     Order 8 (64 coeffs)
    sylph::SHCoefficients visibility(4);
    sylph::SHCoefficients lighting(16);
    sylph::SHCoefficients brdf(8);

    visibility(0, 0) = 0.9;
    lighting(0, 0)   = 2.5;
    brdf(0, 0)       = 0.5;

    // Pack into array of pointers
    std::array<const sylph::SHCoefficients*, 3> inputs = {
        &visibility,
        &lighting,
        &brdf
    };

    sylph::SHGridTransform transform;

    // Compute generalized product truncated to Order 8
    sylph::SHCoefficients output = transform.product(inputs, 8);

    std::cout << "Generalized 3-way product computed. Output order: " << output.order() << "\n";
    return 0;
}
```
