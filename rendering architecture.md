# Sylph: Real-Time Rendering & Graphics Architecture

A comprehensive architectural guide to the computer graphics and real-time rendering systems implemented in Sylph (`include/render/` and `src/render/`).

---

## Table of Contents

1. [High-Level Rendering Architecture](#1-high-level-rendering-architecture)
2. [Precomputed Radiance Transfer (PRT)](#2-precomputed-radiance-transfer-prt)
3. [Object Occlusion Fields (OOF)](#3-object-occlusion-fields-oof)
4. [Shadow Fields: Multiple Dynamic Soft Shadows](#4-shadow-fields-multiple-dynamic-soft-shadows)
5. [Real-Time Glossy Relighting](#5-real-time-glossy-relighting)
6. [Interactive OpenGL Visualizer & SH Lobe Mesh](#6-interactive-opengl-visualizer--sh-lobe-mesh)
7. [System Data Flow & Class Diagram](#7-system-data-flow--class-diagram)

---

## 1. High-Level Rendering Architecture

Sylph is not just an abstract mathematical library; it was engineered to power high-performance real-time computer graphics. 

In traditional rendering, evaluating realistic soft shadows, ambient occlusion, and glossy environment relighting requires casting hundreds of rays per pixel (Monte Carlo path tracing), which is computationally expensive for interactive frame rates.

Sylph solves this by projecting light transport into the Spherical Harmonics basis:
1. **Static scene light transport** is precomputed once onto vertex meshes (**PRT**).
2. **Dynamic moving occluders** carry their own precomputed volumetric shadow fields (**OOF**).
3. **At runtime**, the renderer combines dynamic HDR skies, moving characters, and material BRDFs using **Sylph's microsecond generalized SH product engine**, achieving dynamic soft shadows and glossy reflections at hundreds of frames per second with zero runtime ray tracing!

```
                         THE SYLPH RENDERING PIPELINE
                         
   [ Offline Baking Phase ]                    [ Real-Time Shading Phase ]
   
   Triangle Mesh (.obj)                         Dynamic HDR Environment Light
          |                                                   |
          +---> PRTBaker  ==> PRTMesh (Vertex Transfer)        |
          |                                                   v
          +---> OOFBaker  ==> OOF (3D Volumetric Grid)   ==> ShadowField / GlossyRelighting
                                (Carried by moving obj)       | (Evaluates N-way Product)
                                                              v
                                                       Final Pixel Color
```

---

## 2. Precomputed Radiance Transfer (PRT)

### 2.1 Theoretical Foundation (Sloan et al. 2002)
The rendering equation for diffuse surfaces under distant environment lighting is:
$$L_o(\mathbf{x}) = \frac{\rho}{\pi} \int_{S^2} L_{\text{env}}(\omega) V(\mathbf{x}, \omega) \max(0, \mathbf{n}_{\mathbf{x}} \cdot \omega)\,d\Omega$$
where:
- $L_{\text{env}}(\omega)$ is the incident lighting from the distant sky.
- $V(\mathbf{x}, \omega)$ is the self-visibility (1 if unoccluded, 0 if blocked by other parts of the mesh).
- $\max(0, \mathbf{n} \cdot \omega)$ is the Lambertian cosine factor.

Notice that the product $T(\mathbf{x}, \omega) = V(\mathbf{x}, \omega) \max(0, \mathbf{n}_{\mathbf{x}} \cdot \omega)$ depends **only on the mesh geometry**, not on the lighting!

By precomputing the SH projection of this transfer function at every vertex:
$$\mathbf{t}_v = \int_{S^2} T(\mathbf{x}_v, \omega) \mathbf{y}(\omega)\,d\Omega$$
At runtime, for any dynamic lighting environment $\mathbf{l}_{\text{env}}$, the reflected color at vertex $v$ reduces to an instant **dot product**:
$$L_o(\mathbf{x}_v) = \frac{\rho}{\pi} \left( \mathbf{l}_{\text{env}} \cdot \mathbf{t}_v \right)$$

### 2.2 Classes & Implementation

#### `render::PRTMesh` (`include/render/prt_mesh.hpp`)
Represents a 3D triangle mesh where every vertex carries:
- 3D Position (`Eigen::Vector3f`)
- Normal vector (`Eigen::Vector3f`)
- Precomputed SH transfer coefficients (`sylph::SHCoefficients`, typically Order 4 = 16 coefficients per color channel).

#### `render::PRTBaker` (`include/render/prt_baker.hpp`)
Performs offline CPU Monte Carlo raycasting:
```cpp
static std::shared_ptr<PRTMesh> bake(
    const MeshData& mesh, 
    int order = 4, 
    int num_rays = 1000,
    bool shadows = true
);
```
- For each vertex, casts `num_rays` uniformly distributed over the hemisphere.
- Tests ray-triangle intersections against the mesh to compute self-shadowing $V(\mathbf{x}, \omega)$.
- Accumulates results into vertex SH coefficients.

---

## 3. Object Occlusion Fields (OOF)

### 3.1 The Concept (Zhou et al. 2005 & Paper Section 7.3.2)
While standard PRT handles static self-shadows, it cannot handle **moving, dynamic objects** (like a character walking across a room).

An **Object Occlusion Field (OOF)** solves this by wrapping a 3D bounding box around the dynamic object and discretizing it into a 3D voxel grid.  
**Every single voxel in this grid stores the Spherical Harmonic visibility function of the object as viewed from that point in space.**

```
                              OBJECT OCCLUSION FIELD (OOF)
                              
             +-----------------------------------------------+
             | Voxel (0,2)       Voxel (1,2)     Voxel (2,2) |
             | [ SH: 80% open ]  [ SH: 50% open] [ SH: 90% ] |
             |                                               |
             | Voxel (0,1)       +-------+       Voxel (2,1) |
             | [ SH: 40% open ]  | BUNNY |       [ SH: 40% ] |
             |                   +-------+                   |
             | Voxel (0,0)       Voxel (1,0)     Voxel (2,0) |
             | [ SH: 70% open ]  [ SH: 30% open] [ SH: 70% ] |
             +-----------------------------------------------+
                  Grid follows the Bunny as it moves in 3D!
```

- When the object moves or rotates, the OOF rigid-body transforms with it.
- Any receiver surface outside the object can query the OOF at its world position to determine exactly which angles of the sky the dynamic object is blocking!

### 3.2 Classes & Implementation

#### `render::OOF` (`include/render/oof.hpp`, `src/render/oof.cpp`)
Stores the 3D grid and provides spatial queries:
```cpp
class OOF {
public:
    OOF(int res_x, int res_y, int res_z, 
        const Eigen::Vector3f& min, const Eigen::Vector3f& max, 
        int visibility_order);

    // Queries the visibility function at any continuous 3D point in space
    const Visibility lookup(const Eigen::Vector3f& position) const;
};
```
- **Trilinear SH-Space Interpolation:** In `lookup()`, the 3D position is mapped to fractional voxel coordinates. The method identifies the 8 surrounding voxel corners and performs trilinear interpolation directly on the SH coefficient vectors, eliminating blocky stair-stepped voxel shadow edges and producing smooth dynamic penumbras.

#### `render::OOFBaker` (`include/render/oof_baker.hpp`, `src/render/oof_baker.cpp`)
Bakes the 3D field from a triangle mesh:
```cpp
static std::unique_ptr<OOF> bake(
    const MeshData& mesh,
    int resolution,       // e.g., 16 for a 16x16x16 grid
    int visibility_order, // e.g., 4 or 8
    int num_rays = 500
);
```
- Shoots rays outward in all directions from every voxel center to determine mesh occlusion and projects the result into `sylph::SHCoefficients`.

---

## 4. Shadow Fields: Multiple Dynamic Soft Shadows

### 4.1 Theoretical Foundation (Paper Section 7.3.2)
When multiple dynamic objects (e.g., two moving characters) cast shadows onto a receiver surface, the total received light depends on the **simultaneous product of all occluders**:

$$I(\mathbf{x}) = \int_{S^2} L_{\text{env}}(\omega) \times V_{\text{self}}(\mathbf{x}, \omega) \times V_{\text{OOF}_1}(\mathbf{x}, \omega) \times \dots \times V_{\text{OOF}_k}(\mathbf{x}, \omega)\,d\Omega$$

This is an **$N$-way multiple product integral**.  
Prior methods either had to re-raytrace the scene or recursively multiply and truncate SH coefficients, which caused severe blurring and error accumulation.

### 4.2 Implementation: `render::ShadowField` (`include/render/shadow_field.hpp`)

`ShadowField` evaluates dynamic multi-object soft shadows in real time:
```cpp
Eigen::Vector3f ShadowField::evaluate(
    const Lighting& lighting,
    const Visibility& self_visibility,
    const Eigen::Vector3f& position,
    std::span<const OOF* const> oofs
) const;
```

#### The Execution Flow:
1. Queries each active dynamic occluder's `OOF::lookup(position)` to retrieve its interpolated visibility coefficients at the receiver point.
2. Collects the pointers: `[Lighting, SelfVisibility, OOF_1, OOF_2, ...]`.
3. Calls `sylph::SHGridTransform::product_integral(inputs)`.
4. Sylph's grid transform converts all inputs into a single combined spatial grid, multiplies them pointwise, and evaluates the spherical integral in **microseconds**!

---

## 5. Real-Time Glossy Relighting

### 5.1 Theoretical Foundation (Paper Section 7.3.1)
For non-diffuse surfaces (e.g., polished plastic, ceramic, car paint), the reflected light depends on the view direction $\omega_o$ and the surface BRDF $f_r(\omega_i, \omega_o)$:

$$L_o(\mathbf{x}, \omega_o) = \int_{S^2} L_{\text{env}}(\omega) \times f_r(\omega, \omega_o) \times V(\mathbf{x}, \omega)\,d\Omega$$

This is a **Triple Product Integral** involving three distinct functions:
1. Incident lighting $L_{\text{env}}$ (typically high order, e.g., $N = 16$).
2. Material BRDF $f_r$ (specular lobe centered along reflection vector, order $N = 8$).
3. Visibility $V$ (order $N = 4$).

### 5.2 Implementation: `render::GlossyRelighting` (`include/render/glossy_relighting.hpp`)
```cpp
Eigen::Vector3f GlossyRelighting::evaluate(
    const Lighting& lighting,
    const BRDF& brdf,
    const Visibility& visibility
) const;
```
- Leverages the paper's **generalized product**: inputs with different orders ($n_1 = 16, n_2 = 8, n_3 = 4$) are processed directly without zero-padding overhead.
- Evaluates exact glossy specular reflections with dynamic self-shadowing without relying on prefiltered split-sum approximations.

---

## 6. Interactive OpenGL Visualizer & SH Lobe Mesh

To visualize and debug spherical harmonics in real time, Sylph includes a complete modern OpenGL 4.3 visualization engine.

### 6.1 `render::SHLobeMesh` (`include/render/sh_lobe_mesh.hpp`)
Deforms a high-resolution 3D UV sphere to visually represent any arbitrary `sylph::SHCoefficients`:
- **Geometry Deformation:** For each vertex on the sphere with spherical angles $(\theta, \phi)$, the vertex position is scaled radially by the absolute value of the SH function:
  $$r(\theta, \phi) = |F(\theta, \phi)|$$
- **Color Coding:** 
  - Green / Blue: Positive regions of the harmonic ($F(\theta, \phi) > 0$).
  - Red / Orange: Negative regions of the harmonic ($F(\theta, \phi) < 0$).
- Updates dynamically on the CPU/GPU as coefficients change in real time.

### 6.2 Supporting Engine Infrastructure
- **`render::Renderer` (`renderer.hpp`):** Manages OpenGL state, viewport, depth buffering, and draw calls.
- **`render::Shader` (`shader.hpp`):** Compiles and links GLSL vertex and fragment shaders with automatic uniform caching.
- **`render::Camera` & `CameraController` (`camera.hpp`, `camera_controller.hpp`):** 6-DOF orbital and first-person camera controller handling GLFW keyboard and mouse inputs.
- **`render::ObjLoader` (`obj_loader.hpp`):** Wavefront `.obj` file parser for loading 3D meshes (Stanford Bunny, teapot, sphere).

---

## 7. System Data Flow & Class Diagram

```
+-----------------------------------------------------------------------------------+
|                                 CORE SYLPH MATH                                   |
|                                                                                   |
|   +-------------------+      +--------------------+      +--------------------+   |
|   |  SHCoefficients   |      |   ZHCoefficients   |      |  SHGridTransform   |   |
|   |  Flat N^2 doubles |      |   O(N) Axial Lobe  |      |  Product Engine    |   |
|   +-------------------+      +--------------------+      +--------------------+   |
+-------------|--------------------------|----------------------------|-------------+
              |                          |                            |
              v                          v                            v
+-----------------------------------------------------------------------------------+
|                            RENDERING & GRAPHICS MODULE                            |
|                                                                                   |
|   [ Offline Baking ]                                                              |
|   +--------------+      bakes      +---------------+                              |
|   |   PRTBaker   | --------------> |    PRTMesh    | (Vertex Transfer Vectors)    |
|   +--------------+                 +---------------+                              |
|                                                                                   |
|   +--------------+      bakes      +---------------+                              |
|   |   OOFBaker   | --------------> |      OOF      | (3D Volumetric Grid)         |
|   +--------------+                 +---------------+                              |
|                                            |                                      |
|                                            | queries lookup()                     |
|   [ Real-Time Shading ]                    v                                      |
|   +----------------------+         +---------------+                              |
|   |     ShadowField      | <------ |  Visibility   |                              |
|   | Multi-Product Integr |         +---------------+                              |
|   +----------------------+                                                        |
|              |                                                                    |
|              +--------------------------------------------------+                 |
|                                                                 v                 |
|   +----------------------+                             +------------------+       |
|   |   GlossyRelighting   |                             |     Renderer     |       |
|   | Triple Product Integ |                             |  (OpenGL 4.3)    |       |
|   +----------------------+                             +------------------+       |
|                                                                 ^                 |
|   [ Visualization ]                                             |                 |
|   +----------------------+                                      |                 |
|   |      SHLobeMesh      | -------------------------------------+                 |
|   | (3D Deformed Sphere) |                                                        |
|   +----------------------+                                                        |
+-----------------------------------------------------------------------------------+
```
