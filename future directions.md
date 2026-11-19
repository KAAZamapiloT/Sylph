# Sylph: Future Research Directions and Optimization Roadmap

This document outlines the theoretical, algorithmic, hardware, and rendering research directions extending the **Sylph Generalized Spherical Harmonics library** and the foundational paper:

> **"Generalized Spherical Harmonics Products using Spherical Grids"**  
> *Di An, Jiaqi Wu, Bowen Xu, Lingqi Yan, and Kun Xu*  
> **ACM Transactions on Graphics (TOG 2026), Article 63.**

---

## Roadmap Overview

```
                                      SYLPH FUTURE DIRECTIONS
                                                 |
         +---------------------------------------+---------------------------------------+
         |                                       |                                       |
  [ 1. Core Paper Math ]               [ 2. Hardware & SIMD ]              [ 3. Rendering Systems ]
  - Direct Product Integral (p=1)       - AVX2/FMA Vectorization            - Exact PBR Glossy Relighting
  - Adaptive Band-Truncation (L_act)    - GPU Compute Shaders (GLSL)        - Lazy OOF Dirty Cell Rebaking
  - Mixed-Order Frequency Pruning       - Non-blocking OpenMP Threading     - Multi-Res Sparse Octree OOF
  - Polar Redundancy Mitigation         - Cache-aligned Memory Layouts      - Zero-Copy Trilinear SH Lerp
  - Small-Order Gaunt Tensor (N<=3)
  - Ivanic-Ruedenberg Fast Rotation
```

---

## 1. Core Mathematical Extensions to the Paper

### 1.1 Direct Spatial Quadrature Product Integral & Dot-Product Identity (Paper Section 4 & Section 7.1)

#### The Research Context
In computer graphics (PRT diffuse lighting, ambient occlusion, dynamic shadow fields), the overwhelming majority of SH product operations are **product integrals** where the output order is $p = 1$:
$$I = \int_{S^2} \prod_{k=1}^K F_k(\omega)\,d\omega = \sqrt{4\pi}\,g_{0,0}$$
Currently, running `SHGridTransform::product_integral()` executes the full forward grid transform, multiplies the grid, and then runs the **entire inverse transform** (Eq. 37 azimuthal inverse FFT + Eq. 30 inverse polar Legendre recurrence) only to read coefficient $g_{0,0}$.

#### Formulation: Case $K = 2$ (Analytic Orthonormality)
By SH basis orthonormality ($\int_{S^2} y_l^m y_{l'}^{m'}\,d\Omega = \delta_{l,l'} \delta_{m,m'}$):
$$\int_{S^2} F_1(\omega) F_2(\omega)\,d\Omega = \sum_{l=0}^{\min(n_1, n_2)-1} \sum_{m=-l}^l f^{(1)}_{l,m} f^{(2)}_{l,m} = \mathbf{f}_1 \cdot \mathbf{f}_2$$
- **Complexity:** Drops from $O(N^3)$ grid transform to an $O(N^2)$ flat vector dot product.
- **Expected Runtime:** **$< 10\text{ nanoseconds}$** (compared to $19.7\ \mu\text{s}$ on the grid), a **$2,000\times$ speedup**.

#### Formulation: Case $K \ge 3$ (Direct Spatial Quadrature)
When $p = 1$, only mode $(l = 0, m = 0)$ is required:
1. Azimuthal basis for $m = 0$ is $R_0(\phi) = 1$ (Eq. 21):
   $$s_0(\theta_i) = \frac{2\pi}{N_\phi} \sum_{j=0}^{N_\phi - 1} G(\theta_i, \phi_j)$$
2. Polar basis for $(0, 0)$ is $C_0^0(\theta) = \frac{1}{\sqrt{4\pi}}$ (Eq. 21):
   $$g_{0,0} = \frac{1}{\sqrt{4\pi}} \sum_{i=0}^{N_\theta - 1} w_i s_0(\theta_i)$$

Multiplying by $\sqrt{4\pi}$ yields the exact integral directly from the spatial grid samples:
$$I = \sqrt{4\pi}\,g_{0,0} = \frac{2\pi}{N_\phi} \sum_{i=0}^{N_\theta - 1} w_i \sum_{j=0}^{N_\phi - 1} G(\theta_i, \phi_j)$$
- **Impact:** **100% of the inverse pass is eliminated.** Zero inverse FFTs, zero inverse polar Legendre projections, and zero output coefficient allocations.
- **Expected Speedup:** Instantly cuts $K$-way product integral runtime by **$\approx 50\%$**.

---

### 1.2 Adaptive Band-Truncated / Sparse SH Transforms (Paper Section 6)

#### The Research Context
Real-world spherical functions (environment lighting, visibility fields, diffuse cosine lobes) exhibit rapidly decaying power spectra with increasing degree $l$:
$$\mathcal{P}(l) = \sum_{m=-l}^l |c_{l,m}|^2$$
For smooth signals, bands $l \ge 4$ typically hold less than $0.1\%$ of the total signal energy. However, current representations allocate a flat array of $N^2$ doubles and run all quadrature loops up to nominal order $N$.

#### Mathematical Formulation
Define an energy threshold $\varepsilon$ (e.g., $10^{-6}$) and compute the active band-limit $L_{\text{act}}$:
$$L_{\text{act}} = \max\left\{ l \in [0, N-1] : \frac{\mathcal{P}(l)}{\sum_{k=0}^{N-1} \mathcal{P}(k)} > \varepsilon \right\} + 1$$

1. **Dynamic Quadrature Grid Sizing (Paper Eq. 34):**
   The exact quadrature degree $t$ can be computed using the active orders rather than the nominal orders:
   $$t = \sum_{k=1}^K (L_{\text{act}, k} - 1) + (p - 1) + 1$$
   If two Order-16 inputs have energy concentrated in bands $l < 4$, $t$ drops from $46$ down to $10$.
   - The required Gauss–Legendre grid shrinks from $24 \times 64$ ($1,536$ points) down to $5 \times 16$ ($80$ points).
   - **Throughput Gain:** An immediate **$\sim 15\times$ speedup** with negligible perceptual or numerical error.
2. **Polar Synthesis Early-Exit (Paper Eq. 25):**
   $$s_m(\theta_i) = \sum_{l=|m|}^{L_{\text{act}}-1} f_l^m C_l^m(\theta_i)$$
   Terminates at $L_{\text{act}}$ rather than sweeping all nominal bands.

---

### 1.3 Mixed-Order Azimuthal Frequency Pruning (Paper Section 4 & Eq. 25)

#### The Research Context
In generalized products $\otimes_k^p(\mathbf{f}^{(n_1)}_1, \dots, \mathbf{f}^{(n_k)}_k)$, each input function has a different order $n_i$. For instance, in real-time PRT:
- Visibility: $n_1 = 4$
- Lighting: $n_2 = 16$
- BRDF: $n_3 = 8$

The common grid requires $N_\phi \ge n_G + p - 1 = 26$.

#### Mathematical Formulation
When synthesizing $F_1$ ($n_1 = 4$), the coefficients $f_{l,m}^{(1)}$ are identically zero for all $|m| \ge 4$.
- Out of the $2n_G - 1 = 51$ frequency modes, only $2n_1 - 1 = 7$ modes are non-zero.
- The forward polar loop only needs to evaluate $|m| < n_1$, skipping $86\%$ of the polar synthesis work for the low-order input.
- In the forward IFFT, only the first $n_1$ frequencies are populated; the rest are zero-padded, enabling pruned-spectrum FFT algorithms.

---

### 1.4 Resolving Polar Grid Redundancy (Paper Section 8 "Grid Efficiency")

#### The Research Context
In Section 8 of the paper, the authors explicitly note:
> *"Currently, our method utilizes the SphGrid representation, which requires a grid resolution of $n \times (2n-1)$ to fully reconstruct order-$n$ SH coefficients... The oversampling creates data redundancy, which is most obvious near the poles where points are packed much more densely than at the equator."*

#### Mathematical Formulation: Latitude-Adaptive Bandwidth
Associated Legendre polynomials satisfy:
$$\tilde{P}_l^m(\cos\theta) \propto (\sin\theta)^m$$
At polar nodes $\theta_i$ near the poles ($\theta \to 0$ or $\pi$), $\sin\theta_i \approx 0$. For higher $m$, $\tilde{P}_l^m(\cos\theta_i)$ underflows below machine precision $\varepsilon \approx 10^{-16}$.
- **Cutoff Condition:** For a given latitude node $\theta_i$, prune all azimuthal modes $m$ satisfying:
  $$(\sin\theta_i)^m < \varepsilon_{\text{machine}} \implies m > \frac{\ln \varepsilon_{\text{machine}}}{\ln \sin\theta_i}$$
- **Result:** Polar rings near the poles only evaluate low-$m$ modes ($m \le 2$), eliminating redundant Legendre arithmetic and cutting FFT sizes at high latitudes.

---

### 1.5 Analytical Sparse Gaunt Tensor Crossover for Low Orders ($N \le 3$)

#### The Research Context
The paper replaces Clebsch–Gordan / Gaunt tensor methods because general triple products scale as $O(N^5)$. However, for very low orders, the Gaunt tensor $\Gamma_{ijk} = \int_{S^2} Y_i Y_j Y_k d\Omega$ is extraordinarily sparse:
- **Order 1 ($N = 1$):** 1 coefficient; scalar multiplication.
- **Order 2 ($N = 2$):** 4 coefficients; only 14 non-zero entries in $\Gamma$.
- **Order 3 ($N = 3$):** 9 coefficients; only 114 non-zero entries in $\Gamma$.

#### Implementation Strategy
Implement compile-time unrolled sparse tensor contractions for $N \le 3$:
```cpp
if (max_order <= 3) {
    return sparse_gaunt_product(a, b); // ~25 nanoseconds
} else {
    return product_grid_cached(a, b);   // Microsecond grid transform
}
```
This gives Sylph optimal performance across all orders, outperforming grid setup overhead for low-order graphics workloads.

---

### 1.6 Fast SH Rotation via Direct Recurrence (Ivanic–Ruedenberg / Blanco)

#### The Research Context
In the paper's PRT applications, light sources, surface normals, and BRDF lobes must be rotated dynamically. Rotating an SH function by 3D rotation $\mathbf{R} \in \mathrm{SO}(3)$ transforms coefficients via Wigner D-matrices:
$$c'_{l,m} = \sum_{m'=-l}^l D^l_{m,m'}(\mathbf{R})\,c_{l,m'}$$
Evaluating Wigner D-matrices via Euler angle trigonometry or naive sampling is $O(N^4)$ and numerically unstable at high orders.

#### Formulation
Implement the **Ivanic–Ruedenberg / Blanco recurrence**, which computes the $(2l+1) \times (2l+1)$ matrix $D^l(\mathbf{R})$ directly from $D^{l-1}(\mathbf{R})$ and $D^1(\mathbf{R})$ in $O(N^3)$ total operations:
- Zero trigonometric function calls after band 1.
- Purely algebraic recurrence with zero matrix inversions.
- Provides exact, real-time SH rotation to align lobes before product evaluation.

---

## 2. Hardware Acceleration, SIMD, and GPU Compute

### 2.1 AVX2 / AVX-512 & FMA Vectorization of `PolarBasisCache`

#### The Opportunity
With `PolarBasisCache` already storing normalized Associated Legendre polynomials in contiguous column-major memory, the polar loops in `forward_polar_cached` and `inverse_polar_cached` are structured as:
$$E_{i,m} = \sum_k c_{m+k}^m \cdot \tilde{P}_{m+k}^m(x_i)$$
This is a dense matrix-vector dot product.

#### Implementation
Implement explicit 256-bit AVX2/FMA vector intrinsics:
```cpp
// Process 4 double-precision floats per cycle using FMA
__m256d acc = _mm256_setzero_pd();
for (int k = 0; k < len; k += 4) {
    __m256d c = _mm256_loadu_pd(&coeff[k]);
    __m256d p = _mm256_load_pd(&basis[k * nNorth + i]);
    acc = _mm256_fmadd_pd(c, p, acc);
}
```
- Eliminates compiler autovectorization uncertainty.
- Delivers an immediate **$2\times\text{–}3\times$ speedup** on the polar projection pass on x86_64 CPUs.

---

### 2.2 GPU Compute Shaders (OpenGL 4.3 / CUDA)

#### The Architecture
In scenes with tens of thousands of vertices or massive volumetric grids, evaluating SH transforms sequentially on the CPU becomes bandwidth-limited.
Offloading the transform to GPU Compute Shaders:
1. **Pass 1 (Polar Synthesis Kernel):** Each workgroup handles one latitude ring $\theta_i$. Threads within the workgroup compute $s_m(\theta_i)$ across azimuthal frequencies $m$.
2. **Pass 2 (Azimuthal Synthesis / FFT Kernel):** Compute parallel 1D FFTs across all latitude rings.
3. **Pass 3 (Pointwise Multiply Kernel):** Trivial 2D buffer element-wise multiplication: $G_{i,j} = \prod_k F_{k, i, j}$.
4. **Pass 4 (Inverse Quadrature Kernel):** Sum the product grid into output coefficients $g_l^m$.

On an RTX 4070 (5,888 CUDA cores), this allows evaluating hundreds of thousands of dynamic SH products per frame in parallel.

---

## 3. Computer Graphics & Rendering Applications

### 3.1 Exact PBR Glossy BRDF Relighting (Replacing the Split-Sum Approximation)

#### The Problem
Real-time game engines (Unreal Engine 4/5, Unity HDRP) approximate glossy specular reflection using the Karis 2013 split-sum approximation:
$$\int_\Omega L_i(\omega) f_r(\omega, \omega_o) (\omega \cdot \mathbf{n})\,d\omega \approx \left(\int_\Omega L_i(\omega) D(\omega)\,d\omega\right) \left(\int_\Omega f_r(\omega, \omega_o) (\omega \cdot \mathbf{n})\,d\omega\right)$$
This decouples the lighting from the BRDF, pre-baking the lighting into prefiltered mipmaps and the BRDF into a 2D lookup texture (LUT). While fast, it causes noticeable distortion for anisotropic surfaces, dynamic shadowing, and rough metals.

#### The Sylph Extension
With Sylph's Order-16 SH product running in **$19.7\ \mu\text{s}$** (and symmetric lobes in **$3.5\ \mu\text{s}$**), we can compute the **exact product integral** dynamically:
```cpp
// Per-probe or per-vertex exact reflection
ZHCoefficients brdf_lobe = ggx_to_zh(roughness, view_dir, normal, order);
double specular_radiance = grid_transform.product_integral(env_lighting, brdf_lobe);
```
- **Mathematically exact:** No split-sum decoupling artifacts.
- **Zero texture LUT memory:** Eliminates prefiltered environment mipmaps and 2D BRDF LUTs.

---

### 3.2 Dynamic Lazy OOF Propagation & Spatial Dirty Cell Rebaking

#### The Problem
In the Object Occlusion Field (OOF) volumetric shadow grid, each cell stores an SH visibility function. When an occluder moves, re-raycasting all $R_x \times R_y \times R_z$ voxels every frame requires dozens of milliseconds.

#### The Extension: Bounding Solid-Angle Dirty Tagging
1. When an occluder moves by $\Delta \mathbf{x}$ or rotates by $\Delta \mathbf{q}$, compute the maximum solid-angle change $\Delta \Omega$ for each voxel:
   $$\Delta \Omega \approx \frac{A_{\text{proj}}}{d^2} \|\Delta \mathbf{x}\|$$
2. Only mark cells with $\Delta \Omega > \varepsilon$ as dirty. Far-field cells are unaffected by small motions and skip rebaking.
3. In `ShadowField::evaluate()`, only rebake dirty voxels that are actively queried by receiver geometry in the current camera frustum.

---

### 3.3 Multi-Resolution Sparse Octree OOF

#### The Problem
A uniform $32^3$ OOF grid at Order 16 consumes:
$$32 \times 32 \times 32 \times 256 \times 8\text{ bytes} \approx \mathbf{67.1\text{ MB per occluder}}$$
The vast majority of cells are either completely inside the object (visibility $\equiv 0$) or far away in free space (visibility $\equiv 1$), wasting memory on flat SH coefficient vectors.

#### The Extension
Build a sparse octree that subdivides adaptively based on visibility gradient $\|\nabla V\|$:
- Interior leaf nodes: store a single scalar flag ($0.0$).
- Exterior leaf nodes: store a single scalar flag ($1.0$).
- Boundary / penumbra nodes: store high-order SH coefficients only near silhouette surfaces.
- **Impact:** Memory drops from $O(R^3)$ down to $O(\text{Surface Area})$, shrinking OOF storage from $67\text{ MB}$ to under $5\text{ MB}$ per mesh.

---

### 3.4 Zero-Copy Trilinear SH-Space Interpolation

#### The Optimization
Fix the dangling pointer hazard in `ShadowField::evaluate` when querying multiple dynamic occluders by replacing temporary value copies with a preallocated scratch vector, and vectorize the 8-corner trilinear blend:
$$c_{l,m}(\mathbf{x}) = \sum_{k=0}^7 w_k(\mathbf{x}) c_{l,m}^{(k)}$$
using AVX2 FMA instructions, providing smooth, artifact-free dynamic soft shadow gradients at maximum frame rates.

---

## 4. Prioritized Action Matrix

| Priority | Feature / Research Extension | Primary Benefit | Expected Impact | Code Location |
|---|---|---|---|---|
| **P1** | **Direct Product Integral ($p = 1$)** | Shading & PRT Runtime | **$2,000\times$** ($K=2$), **$2\times$** ($K \ge 3$) | `sh_grid_transform.hpp` |
| **P2** | **Adaptive Band Truncation ($L_{\text{act}}$)** | Dynamic $t$-Sizing & Memory | **$5\times\text{–}15\times$** for smooth lobes | `sh_coefficients.hpp`, `sh_grid_transform.hpp` |
| **P3** | **Mixed-Order Frequency Pruning** | Generalized Products | $1.5\times\text{–}2\times$ forward speedup | `sh_grid_transform.hpp` |
| **P4** | **AVX2 / FMA Vectorization** | Polar Basis Throughput | $2\times\text{–}3\times$ polar speedup | `polar_basis_cache.hpp` |
| **P5** | **Sparse Gaunt Fast-Path ($N \le 3$)** | Small-Order Products | $\sim 30\text{ ns}$ execution | `sh_grid_transform.hpp` |
| **P6** | **Exact PBR Glossy BRDF Reflection** | Shading Quality | Eliminates split-sum error | `glossy_relighting.cpp` |
| **P7** | **Sparse Octree OOF** | Volumetric Memory | $10\times\text{–}15\times$ memory reduction | `oof.hpp`, `oof.cpp` |
| **P8** | **GPU Compute Shaders** | Massive Scene Scaling | Real-time transforms for 100k+ verts | `renderer.cpp`, GLSL kernels |
