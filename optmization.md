# Sylph Optimization Log

A living research and engineering log tracking mathematical, algorithmic, data-structure, and systems optimizations applied to the Sylph Spherical Harmonics library.

---

## Optimization Roadmap & Status

| # | Optimization | Category | Target Complexity / Effect | Status |
|---|--------------|----------|----------------------------|--------|
| 1 | Sparse / Adaptive Truncated SH Coefficients | Data Structure | $O(N^2) \to O(N_{\text{active}}^2)$ storage & loops | Planned |
| 2 | Clenshaw / 3-Term Legendre Recurrence | Algorithm | Eliminate `std::tgamma()` in `Legendre::Pl` (3–8× end-to-end speedup) | Planned (Next Bottleneck) |
| 3 | Lazy OOF Propagation (Dirty Cells) | Algorithm | On-demand partial rebake for small occluder motions | Planned |
| 4 | **FFT Azimuthal Projection & Synthesis** | **Numerical Decomposition** | **$O(N_\theta N_\phi M) \to O(N_\theta N_\phi \log N_\phi)$ (5×–164× azimuthal speedup)** | **Implemented** |
| 5 | GPU Compute Shaders for Forward/Inverse SH | Rendering / GPU | Parallelize $\theta$ rows & $\phi$ samples on GPU | Planned |
| 6 | PBR Split-Sum Replacement via Exact SH Product | Rendering Connection | Exact GGX × Environment integral in SH space | Planned |
| 7 | Trilinear SH-Space OOF Interpolation | Algorithm / Quality | Smooth shadow gradients instead of nearest-voxel steps | Planned |
| 8 | Precomputed Gauss-Legendre Tables | Optimization | Remove Newton iteration from hot path | Planned |
| 9 | Zonal Harmonic (ZH) Fast Path | Algorithm | $O(N^3) \to O(N^2)$ when one input is rotationally symmetric | Planned |
| 10 | Multi-Resolution Sparse Octree OOF | Data Structure | $O(R^3) \to O(S)$ surface-adaptive occlusion grid | Planned |

---

## Log Entry #4 — FFT Azimuthal Projection & Synthesis

- **Date:** 2026-11-18
- **Status:** Implemented & Verified
- **Files Added:**
  - `include/sylph/azimuthal_fft.hpp` — Header-only Radix-2 Cooley-Tukey + Bluestein Chirp-Z FFT with paired-row real-SH projection & synthesis and `thread_local` plan caching.
- **Files Modified:**
  - `include/sylph/sh_grid_transform.hpp` — Replaced direct azimuthal loops in `forward_separable()` (Eq. 20–25) and `inverse_impl()` (Eq. 37) with `AzimuthalFFT`; rounded exact `product()` azimuthal samples to `next_pow2(tt)` for the pure Radix-2 fast path; kept `*_direct` reference variants for apples-to-apples benchmarking.
  - `examples/verification/verification_example.cpp` — Added side-by-side Naive vs. Direct Grid vs. FFT Grid table and an isolated Azimuthal Projection scaling benchmark (orders 4–64, testing both Bluestein and Radix-2).

### 1. Problem & Bottleneck
In the separable SH grid transform:
- **Forward azimuthal synthesis (Eq. 20–25)** evaluated $g(\theta_i, \phi_j) = \sum_{m=-(N-1)}^{N-1} s_{i,m} R_m(\phi_j)$ directly for every $(i, j, m)$.
- **Inverse azimuthal projection (Eq. 37)** evaluated $s_{i,m} = \frac{2\pi}{N_\phi} \sum_{j=0}^{N_\phi - 1} g(\theta_i, \phi_j) R_m(\phi_j)$ directly for every $(i, m, j)$.

Both passes performed $O(N_\theta \cdot N_\phi \cdot M)$ inner-loop `std::cos` / `std::sin` evaluations. Because $\phi_j = \frac{2\pi j}{N_\phi}$ is uniformly spaced around $[0, 2\pi)$, both passes are real-valued Discrete Fourier Transforms (DFT / IDFT) along each $\theta_i$ ring.

### 2. Mathematical Formulation
1. **Inverse Azimuthal Projection $\leftrightarrow$ Forward Complex DFT:**
   For a real row $g_j = g(\theta_i, \phi_j)$, its DFT is:
   $$G_k = \sum_{j=0}^{N_\phi - 1} g_j e^{-i \frac{2\pi k j}{N_\phi}} = \sum_{j=0}^{N_\phi - 1} g_j \cos(k \phi_j) - i \sum_{j=0}^{N_\phi - 1} g_j \sin(k \phi_j)$$
   Matching against the real SH azimuthal basis $R_m(\phi)$ with $r = |m| \bmod N_\phi$:
   $$s_{i,0} = \frac{2\pi}{N_\phi} \operatorname{Re}(G_0), \qquad s_{i,+k} = \frac{2\pi\sqrt{2}}{N_\phi} \operatorname{Re}(G_r), \qquad s_{i,-k} = -\frac{2\pi\sqrt{2}}{N_\phi} \operatorname{Im}(G_r)$$
   The modulo indexing $r = k \bmod N_\phi$ preserves exact aliasing behavior even in low-$t$ approximate mode when $|m| \ge N_\phi / 2$.

2. **Forward Azimuthal Synthesis $\leftrightarrow$ Unnormalized Inverse Complex DFT:**
   Packing the real SH coefficients into a conjugate-symmetric spectrum $S[0 \dots N_\phi - 1]$:
   $$S[0] \mathrel{+}= s_{i,0}, \qquad S[r_+] \mathrel{+}= \frac{s_{i,+k} - i s_{i,-k}}{\sqrt{2}}, \qquad S[r_-] \mathrel{+}= \frac{s_{i,+k} + i s_{i,-k}}{\sqrt{2}}$$
   where $r_+ = k \bmod N_\phi$ and $r_- = (N_\phi - r_+) \bmod N_\phi$, yields $g_j = \sum_{r=0}^{N_\phi - 1} S[r] e^{+i \frac{2\pi r j}{N_\phi}}$ exactly.

3. **Two-Real-Rows-in-One-Complex-FFT Packing:**
   - Projection pairs rows $(2p, 2p+1)$ into $z_j = g^{(2p)}_j + i g^{(2p+1)}_j$, runs a single FFT $Z = \operatorname{DFT}(z)$, and unpacks via conjugate symmetry:
     $$G^{(2p)}_r = \frac{1}{2}(Z_r + \overline{Z_{r'}}), \qquad G^{(2p+1)}_r = \frac{1}{2i}(Z_r - \overline{Z_{r'}})$$
   - Synthesis forms $S = S^{(2p)} + i S^{(2p+1)}$ and reads `real()` / `imag()` directly from one unnormalized IFFT. This halves the total number of FFTs across $\theta$ rows.

4. **Arbitrary $N_\phi$ via Bluestein's Chirp-Z Transform + Power-of-2 Quadrature Padding:**
   - When $N_\phi$ is a power of 2, iterative Radix-2 Cooley-Tukey runs in $O(N_\phi \log N_\phi)$ with precomputed bit-reversal and stage twiddles.
   - When $N_\phi$ is arbitrary (e.g., exact $t = 3N - 2$ or user-supplied approximate $t$), Bluestein's identity $2kj = k^2 + j^2 - (k - j)^2$ converts the length-$N_\phi$ DFT into a power-of-2 cyclic convolution of size $M = \operatorname{next\_pow2}(2N_\phi - 1)$ with precomputed chirp and kernel FFTs.
   - For exact `SHGridTransform::product()`, azimuthal trapezoidal quadrature is exact for *any* $N_\phi \ge t$. Rounding $N_\phi = \operatorname{next\_pow2}(t)$ preserves exactness while always activating the pure Radix-2 fast path (~4× faster than Bluestein).

### 3. Verification & Benchmarks (MSVC x64 `/O2`, i7-13700HX)

- **System-Level Mathematical Validation (`tests/SysVerificationMath.hpp`):** **19 / 19 PASSED** (0 failures), covering orthonormality, round-trip, exact generalized products, product integrals, exact-$t$ threshold, lower-$t$ approximate mode, commutativity, Parseval energy identity, and Eq. 34 boundaries.

#### Isolated Azimuthal Projection (Eq. 37) — Direct vs. FFT Scaling

| Order $N$ | $N_\phi$ (exact) | $N_\phi$ (pow2) | Direct (ms) | FFT Bluestein (ms) | FFT Radix-2 (ms) | Speedup (Radix-2) | Speedup (Bluestein) | Max Error |
|----------:|-----------------:|----------------:|------------:|-------------------:|-----------------:|------------------:|--------------------:|----------:|
| 4 | 10 | 16 | 0.0014 | 0.0011 | 0.0003 | **5.44×** | 1.25× | $2.22 \times 10^{-15}$ |
| 8 | 22 | 32 | 0.0180 | 0.0045 | 0.0012 | **15.50×** | 3.97× | $4.86 \times 10^{-15}$ |
| 12 | 34 | 64 | 0.0744 | 0.0146 | 0.0033 | **22.51×** | 5.09× | $5.55 \times 10^{-15}$ |
| 16 | 46 | 64 | 0.2193 | 0.0210 | 0.0051 | **43.11×** | 10.44× | $9.16 \times 10^{-15}$ |
| 24 | 70 | 128 | 0.7232 | 0.0619 | 0.0141 | **51.18×** | 11.69× | $7.49 \times 10^{-15}$ |
| 32 | 94 | 128 | 1.7635 | 0.0854 | 0.0205 | **85.84×** | 20.64× | $1.02 \times 10^{-14}$ |
| 48 | 142 | 256 | 6.3679 | 0.2852 | 0.0701 | **90.89×** | 22.33× | $1.83 \times 10^{-14}$ |
| 64 | 190 | 256 | 14.3453 | 0.3604 | 0.0874 | **164.17×** | 39.80× | $1.93 \times 10^{-14}$ |

#### End-to-End SH Product — Naive vs. Direct Grid vs. FFT Grid

| Order $N$ | Naive $O(N^5)$ (ms) | Direct Grid (ms) | FFT Grid (ms) | Speedup vs. Naive | FFT / Direct | Max Error (FFT vs. Direct) |
|----------:|--------------------:|-----------------:|--------------:|------------------:|-------------:|---------------------------:|
| 4 | 0.2173 | 0.0331 | 0.0235 | 9.25× | **1.41×** | $2.78 \times 10^{-16}$ |
| 6 | 1.6720 | 0.1370 | 0.1023 | 16.34× | **1.34×** | $1.11 \times 10^{-15}$ |
| 8 | 7.4090 | 0.5369 | 0.3395 | 21.82× | **1.58×** | $2.83 \times 10^{-15}$ |
| 10 | 25.4388 | 1.3114 | 0.9144 | 27.82× | **1.43×** | $4.55 \times 10^{-15}$ |
| 12 | 73.2435 | 2.0993 | 1.7928 | 40.85× | **1.17×** | $5.96 \times 10^{-15}$ |
| 16 | 335.3346 | 8.1110 | 7.1125 | **47.15×** | **1.14×** | $1.44 \times 10^{-14}$ |

### 4. Profiling Takeaway (Amdahl's Law)
- The FFT azimuthal pass is **43× faster at Order 16** and **164× faster at Order 64**, with **no crossover penalty** (already 5.4× faster at Order 4).
- End-to-end product speedup at Order 16 is 1.14× because the remaining **~87% of runtime** is now spent in the polar Legendre pass (`Legendre::Pl` calling `std::tgamma`). Implementing **#2 (3-term Legendre recurrence)** next will unlock the full combined speedup.
