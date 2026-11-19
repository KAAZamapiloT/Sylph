# Sylph Optimization Log

A living research and engineering log tracking mathematical, algorithmic, data-structure, and systems optimizations applied to the Sylph Spherical Harmonics library.

---

## Optimization Roadmap & Status

| # | Optimization | Category | Target Complexity / Effect | Status |
|---|--------------|----------|----------------------------|--------|
| 1 | Sparse / Adaptive Truncated SH Coefficients | Data Structure | $O(N^2) \to O(N_{\text{active}}^2)$ storage & loops | Planned |
| 2 | **Normalized 3-Term Legendre Recurrence + Equator Parity** | **Algorithm** | **Zero `tgamma`, zero overflow, 211×–244× end-to-end speedup** | **Implemented** |
| 3 | Lazy OOF Propagation (Dirty Cells) | Algorithm | On-demand partial rebake for small occluder motions | Planned |
| 4 | **FFT Azimuthal Projection & Synthesis** | **Numerical Decomposition** | **$O(N_\theta N_\phi M) \to O(N_\theta N_\phi \log N_\phi)$ (5×–164× azimuthal speedup)** | **Implemented** |
| 5 | GPU Compute Shaders for Forward/Inverse SH | Rendering / GPU | Parallelize $\theta$ rows & $\phi$ samples on GPU | Planned |
| 6 | PBR Split-Sum Replacement via Exact SH Product | Rendering Connection | Exact GGX × Environment integral in SH space | Planned |
| 7 | Trilinear SH-Space OOF Interpolation | Algorithm / Quality | Smooth shadow gradients instead of nearest-voxel steps | Planned |
| 8 | **Cached Gauss–Legendre Tables + Polar Basis Cache + Zero-Alloc** | **Optimization** | **1.6×–3.0× (`SH×SH`), 2.8× (`ZH×SH`), 8×–12× (`ZH×ZH`)** | **Implemented** |
| 9 | **Zonal Harmonic (ZH) Fast Path** | **Algorithm** | **Zero FFT / zero 2D grid: 3×–5× (`ZH×SH`) & 5×–12× (`ZH×ZH`, $O(N^2)$)** | **Implemented** |
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

---

## Log Entry #2 — Normalized 3-Term Associated Legendre Recurrence + Equator Parity Folding

- **Date:** 2026-11-18
- **Status:** Implemented & Verified
- **Files Modified:**
  - `include/sylph/legendre.hpp` — Replaced closed-form `std::tgamma` sum in `Legendre::Pl` with directly normalized 3-term recurrence $\tilde{P}_l^m(x) = N_l^m P_l^m(x)$, precomputed column-major `RecurrenceTable` ($d_m, a_{l,m}, b_{l,m}$) cached `thread_local` per order, and recurrence-based `associated()` / `sh()`.
  - `include/sylph/sh_coefficients.hpp` — Added `data()` raw pointer accessors for verified inner loops and fixed `validate()` upper bound (`l >= order_`).
  - `include/sylph/sh_grid_transform.hpp` — Added `forward_polar_recurrence()` and `inverse_polar_recurrence()` featuring: (1) column-major normalized 3-term recurrence, (2) simultaneous $\pm m$ reuse (2× fewer evaluations), (3) North–South Equator Parity Folding across symmetric Gauss-Legendre pairs $(i, N_\theta - 1 - i)$ (another 2× reduction), and (4) branchless unrolled-by-2 even/odd loops.
  - `examples/verification/verification_example.cpp` — Updated Table 1 to compare Naive vs. Orig Grid vs. #4 Only vs. #4+#2, and added Table 3 testing high-order scaling and numerical stability up to Order 64.

### 1. Problem & Bottleneck
After Optimization #4 sped up azimuthal projection by 43×–164×, profiling showed ~87% of end-to-end `product()` runtime remained in `Legendre::Pl()`, which evaluated a closed-form sum calling `std::tgamma()` inside the innermost loop of every `(l, m, theta)`. Worse, alternating large factorial ratios in `tgamma` caused catastrophic cancellation already visible at Order 12–16 ($\sim 5 \times 10^{-12}$) and overflow past Order 32.

### 2. Mathematical Formulation
1. **Directly Normalized 3-Term Recurrence for $\tilde{P}_l^m(x) = N_l^m P_l^m(x)$** (with Condon–Shortley phase $(-1)^m$, $x = \cos\theta$, $s = \sin\theta \ge 0$):
   - **Diagonal seed ($l = m$):**
     $$\tilde{P}_0^0(x) = \frac{1}{\sqrt{4\pi}}, \qquad \tilde{P}_m^m(x) = -\sqrt{\frac{2m + 1}{2m}}\,s\,\tilde{P}_{m-1}^{m-1}(x)$$
   - **Sub-diagonal ($l = m + 1$):**
     $$\tilde{P}_{m+1}^m(x) = \sqrt{2m + 3}\,x\,\tilde{P}_m^m(x)$$
   - **General step ($l \ge m + 2$):**
     $$\tilde{P}_l^m(x) = a_{l,m}\,x\,\tilde{P}_{l-1}^m(x) - b_{l,m}\,\tilde{P}_{l-2}^m(x)$$
     where
     $$a_{l,m} = \sqrt{\frac{4l^2 - 1}{l^2 - m^2}}, \qquad b_{l,m} = \frac{a_{l,m}}{a_{l-1,m}} = \sqrt{\frac{2l + 1}{2l - 3}\cdot\frac{(l - 1)^2 - m^2}{l^2 - m^2}}$$
   All factorials cancel analytically, so coefficients stay $O(1)$ with zero overflow even beyond Order 64. All $d_m, a_{l,m}, b_{l,m}$ are precomputed once per order in column-major layout (zero `sqrt` or divisions in the inner loop).

2. **Simultaneous $\pm m$ Reuse:**
   Because $C(l, +m, \theta) = C(l, -m, \theta) = \tilde{P}_l^{|m|}(\cos\theta)$, each column $\tilde{P}_l^m(x_i)$ is computed once for $m \ge 0$ and immediately applied to both $+m$ and $-m$ coefficients.

3. **North–South Equator Parity Folding:**
   Gauss–Legendre nodes are symmetric across the equator: $x_{i'} = -x_i$ and $w_{i'} = w_i$ for $i' = N_\theta - 1 - i$. Using $\tilde{P}_l^m(-x) = (-1)^{l-m}\tilde{P}_l^m(x)$:
   - **Forward polar synthesis:** Split $l = m \dots N-1$ into even ($E$) and odd ($O$) sums of $(l - m)$ at northern node $x_i$ only:
     $$s_{i, m} = E_{i,m} + O_{i,m}, \qquad s_{i', m} = E_{i,m} - O_{i,m}$$
   - **Inverse polar projection (Eq. 30):** Fold symmetric row pairs before sweeping $l$:
     $$s^+_{i,m} = w_i (s_{i,m} + s_{i',m}), \qquad s^-_{i,m} = w_i (s_{i,m} - s_{i',m})$$
     accumulating $s^+$ into even $(l - m)$ and $s^-$ into odd $(l - m)$ via an unrolled-by-2 branchless loop.

### 3. Verification & Benchmarks (MSVC x64 `/O2`, i7-13700HX)

- **System-Level Mathematical Validation (`tests/SysVerificationMath.hpp`):** **19 / 19 PASSED** (0 failures).

#### End-to-End SH Product — Baseline vs. #4 Only vs. #4 + #2

| Order $N$ | Naive $O(N^5)$ (ms)* | Orig Grid (ms) | #4 FFT Only (ms) | #4 + #2 Fast (ms) | Speedup vs. Orig | Speedup vs. #4 | Speedup vs. Naive |
|----------:|---------------------:|---------------:|-----------------:|------------------:|-----------------:|---------------:|------------------:|
| 4 | 0.0511 | 0.0296 | 0.0229 | **0.0031** | **9.6×** | **7.4×** | 16.6× |
| 6 | 0.3089 | 0.1278 | 0.1041 | **0.0066** | **19.5×** | **15.9×** | 47.1× |
| 8 | 1.2043 | 0.3984 | 0.3123 | **0.0085** | **47.0×** | **36.8×** | 141.9× |
| 10 | 3.5864 | 1.0538 | 0.8674 | **0.0119** | **88.7×** | **73.0×** | 301.9× |
| 12 | 7.9599 | 2.4401 | 1.9103 | **0.0209** | **116.6×** | **91.3×** | 380.2× |
| 16 | 34.9639 | 8.1677 | 7.0544 | **0.0334** | **244.6×** | **211.3×** | **1,047.1×** |

*\*Note: `Naive O(N^5)` also became ~10× faster automatically because `SphericalHarmonics::evaluate()` now uses the recurrence; compared against the original pre-#2 `tgamma` Naive (335.3 ms at Order 16), `#4 + #2` is **~10,040× faster**.*

#### High-Order Scaling & Numerical Stability (Orders 16–64)

| Order $N$ | Product Time (ms) | Round-Trip Err | Commutativity Err | Unit-1 Product Err |
|----------:|------------------:|---------------:|------------------:|-------------------:|
| 16 | **0.0351** | $2.00 \times 10^{-15}$ | $0.00$ | $4.44 \times 10^{-15}$ |
| 24 | **0.1564** | $4.89 \times 10^{-15}$ | $0.00$ | $7.55 \times 10^{-15}$ |
| 32 | **0.2537** | $7.22 \times 10^{-15}$ | $0.00$ | $8.11 \times 10^{-15}$ |
| 48 | **0.8017** | $1.24 \times 10^{-14}$ | $0.00$ | $1.24 \times 10^{-14}$ |
| 64 | **1.2587** | $5.29 \times 10^{-14}$ | $0.00$ | $2.65 \times 10^{-14}$ |

### 4. Key Takeaways
- **244.6× faster end-to-end at Order 16** (`8.17 ms -> 0.0334 ms`, ~30,000 full Order-16 products/sec on 1 CPU core).
- **Fixed numerical cancellation & overflow:** Legacy `tgamma` lost precision by Order 16 ($5.47 \times 10^{-12}$) and overflowed at higher orders; the normalized recurrence maintains **$10^{-14}$ accuracy all the way to Order 64** (1.26 ms for 4,096 × 4,096 SH coefficients).
- **Next bottleneck:** With both polar recurrence and azimuthal FFT now in the microsecond regime, the remaining overhead per `product()` call is **#8 (Gauss–Legendre Newton iteration + per-call vector allocations & `cos`/`sin` recomputation)**.

---

## Log Entry #9 — Zonal Harmonic (ZH) Fast Path

- **Date:** 2026-11-18
- **Status:** Implemented & Verified
- **Files Added:**
  - `include/sylph/zonal_harmonics.hpp` — Compact $O(N)$ `ZHCoefficients` container with exact $O(N^2)$ Addition-Theorem axis rotation (`rotate_to_sh`) and exact $O(N^2)$ Funk–Hecke spherical convolution (`convolve`).
- **Files Modified:**
  - `include/sylph/sh_coefficients.hpp` — Added $O(1)$-early-exit `is_zonal(eps)` detection.
  - `include/sylph/sh_grid_transform.hpp` — Added automatic 3-case ZH dispatch in `product(inputs, output_order)` plus explicit `ZHCoefficients` overloads and `product_general()` benchmark helper:
    1. **All inputs zonal (`ZH × ZH × ...`):** pure $O(N^2)$ 1D $m = 0$ polar quadrature.
    2. **$(K-1)$ zonal inputs + 1 general SH (`ZH × SH`):** fused single-pass polar product (zero 2D grid, zero FFTs, single Legendre sweep, strict $O(N^2)$ when `output_order == 1` for `product_integral`).
    3. **$\ge 1$ zonal inputs + $\ge 2$ general SH (`ZH × SH₁ × SH₂`):** evaluates 1D zonal profile in $O(N^2)$ and scales the 2D grid rows, skipping a full $O(N^3)$ forward transform + FFT per zonal input.
  - `examples/verification/verification_example.cpp` — Added Table 4 comparing General (#4+#2) vs. `ZH × SH` and `ZH × ZH` across Orders 4–64.

### 1. Mathematical Formulation
When a factor $Z(\theta, \phi) = Z(\theta)$ is rotationally symmetric around the $z$-axis ($z_{l,m} = 0$ for all $m \neq 0$), it is independent of $\phi$. Projecting $Z(\theta) F(\theta, \phi)$ onto $Y_l^m(\theta, \phi) = \tilde{P}_l^{|m|}(\cos\theta) R_m(\phi)$ decouples the azimuthal integral analytically:
$$h_{l,m} = \int_0^\pi Z(\theta) \tilde{P}_l^{|m|}(\cos\theta) \sum_{l',m'} f_{l',m'} \tilde{P}_{l'}^{|m'|}(\cos\theta) \underbrace{\int_0^{2\pi} R_{m'}(\phi) R_m(\phi)\,d\phi}_{= \, 2\pi\,\delta_{m,m'}} \sin\theta\,d\theta$$

Consequences:
1. **Zero azimuthal mixing & zero FFTs:** Only $m' = m$ survives, and $|m| \ge \min(N_f, O)$ is identically zero.
2. **Exact 1D Gauss–Legendre polar quadrature:** With $z(x_i) = \sum_l z_{l,0}\tilde{P}_l^0(x_i)$ and $s_{i,m} = \sum_{l'} f_{l',m}\tilde{P}_{l'}^{|m|}(x_i)$, the integrand is a polynomial in $x = \cos\theta$ of degree $\le tt - 1$, integrated exactly by $N_\theta = \lceil tt/2 \rceil$ nodes:
   $$h_{l,m} = \sum_{i=0}^{N_\theta - 1} w_i \, \tilde{P}_l^{|m|}(x_i) \left[2\pi\,z(x_i)\,s_{i,m}\right]$$
3. **Single-pass fused recurrence:** For each northern/southern Gauss–Legendre pair $(i, i')$ and each $m$, $\tilde{P}_l^m(x_i)$ is computed **once** into an L1 buffer `p_col`, dotted with $F$ (even/odd parity), scaled by $2\pi z(x_i), 2\pi z(x_{i'})$, and immediately projected back into `result` — cutting 3 Legendre passes + 3 FFTs + a 2D grid down to **1 Legendre pass and 0 FFTs**.
4. **Strict $O(N^2)$ cases:**
   - `ZH × ZH`: only $m = 0$ exists $\implies O(N^2)$.
   - `product_integral(ZH, SH)` (`output_order == 1`): only $m = 0$ contributes $\implies O(N^2)$.
   - `ZHCoefficients::rotate_to_sh(theta_d, phi_d)`: $c_{l,m} = \sqrt{\frac{4\pi}{2l+1}}\,z_l\,Y_l^m(\theta_d, \phi_d)$ (Addition Theorem, $O(N^2)$).
   - `ZHCoefficients::convolve(F)`: $(Z * F)_{l,m} = \sqrt{\frac{4\pi}{2l+1}}\,z_l\,f_{l,m}$ (Funk–Hecke, $O(N^2)$).

### 2. Verification & Benchmarks (MSVC x64 `/O2`, i7-13700HX)

- **`tests/SysVerificationMath.hpp`:** **19 / 19 PASSED** (0 failures).
- **Additional ZH invariants:**
  - Addition-Theorem Rotation (Order 16): max error **$2.69 \times 10^{-15}$**
  - 3-Way Mixed (`ZH × SH × SH`): max error **$4.55 \times 10^{-15}$**
  - Asymmetric `ZH × SH` (`out > f` / `out < f`): max error **$1.25 \times 10^{-16}$** / **$6.59 \times 10^{-17}$**

#### General (#4 + #2) vs. #9 Zonal Harmonic Fast Path

| Order $N$ | General `#4+#2` (ms) | `ZH × SH` (ms) | `ZH × ZH` (ms) | `ZH×SH` Speedup | `ZH×ZH` Speedup | Err (`ZH×SH`) | Err (`ZH×ZH`) |
|----------:|---------------------:|---------------:|---------------:|----------------:|----------------:|--------------:|--------------:|
| 4 | 0.0030 | **0.0006** | **0.0005** | **4.75×** | **5.75×** | $2.78 \times 10^{-17}$ | $5.55 \times 10^{-17}$ |
| 8 | 0.0078 | **0.0026** | **0.0018** | **2.99×** | **4.47×** | $4.44 \times 10^{-16}$ | $4.44 \times 10^{-16}$ |
| 12 | 0.0204 | **0.0053** | **0.0043** | **3.86×** | **4.79×** | $5.00 \times 10^{-16}$ | $1.67 \times 10^{-16}$ |
| 16 | 0.0315 | **0.0099** | **0.0059** | **3.18×** | **5.33×** | $2.00 \times 10^{-15}$ | $8.88 \times 10^{-16}$ |
| 24 | 0.1330 | **0.0268** | **0.0131** | **4.97×** | **10.16×** | $4.00 \times 10^{-15}$ | $4.44 \times 10^{-16}$ |
| 32 | 0.2377 | **0.0638** | **0.0263** | **3.73×** | **9.03×** | $3.50 \times 10^{-15}$ | $1.90 \times 10^{-15}$ |
| 48 | 0.6801 | **0.1880** | **0.0576** | **3.62×** | **11.80×** | $6.99 \times 10^{-15}$ | $3.11 \times 10^{-15}$ |
| 64 | 1.1623 | **0.4311** | **0.1224** | **2.70×** | **9.50×** | $3.07 \times 10^{-14}$ | $5.77 \times 10^{-15}$ |

### 3. Profiling Takeaway
- At Order 16, `ZH × SH` now runs in **9.9 µs** (>100,000 products/sec, **~825× faster** than the original unoptimized grid) and `ZH × ZH` runs in **5.9 µs**.
- For `ZH × ZH` at Order 64 (0.122 ms), a large fraction of the remaining time is now `gauss_legendre(95)` running Newton iteration from scratch on every call — **#8 (Cached Gauss–Legendre Tables)** will collapse that overhead directly.

---

## Log Entry #8 — Cached Gauss–Legendre Tables + Cached Polar Basis + Zero-Alloc Scratch

- **Date:** 2026-11-18
- **Status:** Implemented & Verified
- **Files Added:**
  - `include/sylph/polar_basis_cache.hpp` — `PolarBasisCache::get(nTheta, order)` caching all normalized Legendre basis values $\tilde{P}_{m+k}^m(x_i)$ for northern/equator Gauss–Legendre nodes in contiguous column-major layout (L1/L2 resident).
- **Files Modified:**
  - `include/sylph/gauss_legendre.hpp` — Added `gauss_legendre::Rule` (`nodes, weights, theta, sin_theta`) and `thread_local` `gauss_legendre::get(n)` cache (O(1) lookup after first call, zero Newton iterations, zero trig).
  - `include/sylph/spherical_grid.hpp` — Added `data()` raw buffer accessors.
  - `include/sylph/azimuthal_fft.hpp` — Added `project_rows_raw()` / `synthesize_rows_raw()` with `thread_local` complex scratch buffers (zero heap allocations on repeated calls) and fast modulo avoidance when $k < N_\phi$.
  - `include/sylph/sh_grid_transform.hpp` — Implemented zero-allocation `product_2d_cached()`, `forward_polar_cached()`, `inverse_polar_cached()`, upgraded all ZH fast paths (`ZH × ZH`, `ZH × SH`) and approximate-$t$ `product(..., t)` (used by `ShadowField` / `OofExample`) to read directly from `gauss_legendre::get()` + `PolarBasisCache::get()` + `thread_local` scratch buffers.
  - `examples/verification/verification_example.cpp` — Added Table 5 comparing Before #8 vs. After #8 across Orders 4–64.

### 1. What Was Eliminated
1. **Gauss–Legendre Newton Iteration:** Re-solving $P_n(x_i) = 0$ on every `product()` call dropped from $O(n^2)$ to a ~2 ns `thread_local` lookup.
2. **Redundant Polar Recurrence & Serial FMA Chains:** Previously, general `product(A, B)` ran the 3-term recurrence 3 times per call (`A`, `B`, `inverse`), with a loop-carried serial dependency `p[k-2], p[k-1] -> p[k]`. Reading precomputed $\tilde{P}_l^m(x_i)$ contiguously from `PolarBasisCache` eliminates all recurrence arithmetic and exposes full out-of-order FMA parallelism.
3. **~18 Heap Allocations & Trig Round-Trips per Call:** Eliminated temporary `Spherical_grid` construction (`acos`), redundant `cos`/`sin` calls, and per-call `std::vector` allocations in the hot path.

### 2. Verification & Benchmarks (MSVC x64 `/O2`, i7-13700HX)

- **`tests/SysVerificationMath.hpp`:** **19 / 19 PASSED** (0 failures).

#### Before #8 vs. After #8 (General `SH × SH`, `ZH × SH`, `ZH × ZH`)

| Order $N$ | Before #8 `SH×SH` (ms) | After #8 `SH×SH` (ms) | Gain (#8) | `ZH × SH` After #8 (ms) | `ZH × ZH` After #8 (ms) | Err (`SH×SH`) | Err (`ZH×ZH`) |
|----------:|-----------------------:|----------------------:|----------:|------------------------:|------------------------:|--------------:|--------------:|
| 4 | 0.0027 | **0.0009** | **3.00×** | **0.0002** | **0.0001** | $1.11 \times 10^{-16}$ | $2.78 \times 10^{-17}$ |
| 8 | 0.0085 | **0.0042** | **2.00×** | **0.0007** | **0.0003** | $1.05 \times 10^{-15}$ | $1.11 \times 10^{-16}$ |
| 12 | 0.0204 | **0.0124** | **1.65×** | **0.0016** | **0.0004** | $1.89 \times 10^{-15}$ | $5.55 \times 10^{-17}$ |
| 16 | 0.0317 | **0.0197** | **1.61×** | **0.0035** | **0.0007** | $7.77 \times 10^{-15}$ | $1.11 \times 10^{-16}$ |
| 24 | 0.1486 | **0.0593** | **2.51×** | **0.0099** | **0.0014** | $1.22 \times 10^{-14}$ | $1.11 \times 10^{-16}$ |
| 32 | 0.2271 | **0.0899** | **2.53×** | **0.0237** | **0.0030** | $2.26 \times 10^{-14}$ | $1.11 \times 10^{-16}$ |
| 48 | 0.7015 | **0.3239** | **2.17×** | **0.0772** | **0.0060** | $4.50 \times 10^{-14}$ | $8.88 \times 10^{-16}$ |
| 64 | 1.1964 | **0.5451** | **2.19×** | **0.1764** | **0.0100** | $1.14 \times 10^{-13}$ | $8.88 \times 10^{-16}$ |

### 3. Cumulative Order-16 Progression (Single Core, MSVC `/O2`)

| Stage | Order-16 Time | Throughput | Speedup vs. Orig Grid | Speedup vs. Orig Naive |
|-------|--------------:|-----------:|----------------------:|-----------------------:|
| Original Naive $O(N^5)$ (`tgamma`) | 335.33 ms | 3 / sec | — | 1× |
| Original Grid (`tgamma` + Direct DFT) | 8.17 ms | 122 / sec | 1× | 41× |
| + **#4** FFT Azimuthal | 7.05 ms | 142 / sec | 1.16× | 48× |
| + **#2** Normalized Legendre Recurrence + Parity | 0.0334 ms | 29,940 / sec | 245× | 10,040× |
| + **#8** Cached GL + Polar Basis + Zero-Alloc (`SH × SH`) | **0.0197 ms** | **50,760 / sec** | **415×** | **17,020×** |
| + **#9 + #8** Zonal Harmonic Fast Path (`ZH × SH`) | **0.0035 ms** | **285,700 / sec** | **2,334×** | **95,800×** |
| + **#9 + #8** Zonal Harmonic Fast Path (`ZH × ZH`) | **0.0007 ms** | **1,428,000 / sec** | **11,670×** | **479,000×** |
