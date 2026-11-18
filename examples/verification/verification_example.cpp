#include "verification_example.hpp"

#include "sylph/sh_coefficients.hpp"
#include "sylph/sh_grid_transform.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/gauss_legendre.hpp"
#include "sylph/azimuthal_fft.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <cmath>
#include <vector>
#include <types/constants.hpp>

namespace examples
{

static sylph::SHCoefficients naive_product(
    const sylph::SHCoefficients& a,
    const sylph::SHCoefficients& b,
    int output_order)
{
    int full_order = a.order() + b.order() - 1;
    int tt = full_order + output_order - 1;
    int nTheta = (tt + 1) / 2;
    int nPhi = tt;

    sylph::gauss_legendre gl(nTheta);
    std::vector<double> weights = gl.weights();
    std::vector<double> nodes = gl.nodes();

    sylph::SHCoefficients result(output_order);
    sylph::SphericalHarmonics sh;

    const double dPhi = 2.0 * sylph::pi / nPhi;

    for (int i = 0; i < nTheta; ++i) {
        double theta = std::acos(nodes[i]);
        for (int j = 0; j < nPhi; ++j) {
            double phi = j * dPhi;
            double valA = sh.reconstruct(a, theta, phi);
            double valB = sh.reconstruct(b, theta, phi);

            for (int l = 0; l < output_order; ++l) {
                for (int m = -l; m <= l; ++m) {
                    double y = sh.evaluate(l, m, theta, phi);
                    result(l, m) += valA * valB * y * weights[i] * dPhi;
                }
            }
        }
    }
    return result;
}

static void direct_azimuthal_project(
    const sylph::Spherical_grid& grid,
    int output_order,
    std::vector<double>& s)
{
    const int nTheta = grid.theta_size();
    const int nPhi   = grid.phi_size();
    const int mCount = 2 * output_order - 1;
    const int mZero  = output_order - 1;

    s.assign(static_cast<std::size_t>(nTheta) * static_cast<std::size_t>(mCount), 0.0);
    const double phiScale = 2.0 * sylph::pi / static_cast<double>(nPhi);
    const double sqrt2 = std::numbers::sqrt2_v<double>;

    for (int i = 0; i < nTheta; ++i) {
        for (int m = -output_order + 1; m <= output_order - 1; ++m) {
            double sum = 0.0;
            for (int j = 0; j < nPhi; ++j) {
                const double phi = grid.phi(j);
                double r = 1.0;
                if (m > 0) {
                    r = sqrt2 * std::cos(m * phi);
                } else if (m < 0) {
                    r = sqrt2 * std::sin(std::abs(m) * phi);
                }
                sum += grid(i, j) * r;
            }
            s[static_cast<std::size_t>(i) * static_cast<std::size_t>(mCount)
              + static_cast<std::size_t>(m + mZero)] = phiScale * sum;
        }
    }
}

VerificationExample::VerificationExample(int width, int height)
{
    resize(width, height);
}

void VerificationExample::process_event(const SDL_Event&, SDL_Window*)
{
}

void VerificationExample::update(float)
{
    if (!run_once_) {
        run_verification();
        run_once_ = true;
    }
}

void VerificationExample::resize(int width, int height)
{
    if (height <= 0) return;
    camera_.set_aspect_ratio(static_cast<float>(width) / static_cast<float>(height));
}

void VerificationExample::render(render::Renderer&)
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void VerificationExample::run_verification()
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    sylph::SHGridTransform grid_transform;

    // ========================================================================
    // Table 1: End-to-End SH Product (Naive vs Direct Grid vs FFT Grid)
    // ========================================================================
    std::cout << "\n========================================================================================================\n";
    std::cout << "             SH Product Mathematical Verification & Speedup Analysis (End-to-End)\n";
    std::cout << "========================================================================================================\n";
    std::cout << "Comparing O(N^5) Naive vs O(N^3) Direct Grid vs O(N^3) FFT Azimuthal Grid (Radix-2)\n\n";

    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(14) << "Naive Time"
              << std::setw(14) << "Direct Grid"
              << std::setw(14) << "FFT Grid"
              << std::setw(13) << "vs Naive"
              << std::setw(13) << "FFT/Direct"
              << std::setw(15) << "Err(FFT-Dir)"
              << std::setw(15) << "Err(FFT-Nai)"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------------------------------\n" << std::flush;

    std::vector<int> orders = {4, 6, 8, 10, 12, 16};

    for (int order : orders) {
        sylph::SHCoefficients a(order);
        sylph::SHCoefficients b(order);
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                a(l, m) = dist(gen);
                b(l, m) = dist(gen);
            }
        }

        // Warm up FFT plan cache for this order
        (void)grid_transform.product(a, b, order);

        auto t0 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_naive = naive_product(a, b, order);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_naive = t1 - t0;

        auto t2 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_direct = grid_transform.product_direct(a, b, order);
        auto t3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_direct = t3 - t2;

        auto t4 = std::chrono::high_resolution_clock::now();
        sylph::SHCoefficients res_fft = grid_transform.product(a, b, order);
        auto t5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_fft = t5 - t4;

        double max_err_dir = 0.0;
        double max_err_nai = 0.0;
        for (int l = 0; l < order; ++l) {
            for (int m = -l; m <= l; ++m) {
                max_err_dir = std::max(max_err_dir, std::abs(res_direct(l, m) - res_fft(l, m)));
                max_err_nai = std::max(max_err_nai, std::abs(res_naive(l, m) - res_fft(l, m)));
            }
        }

        const double speedup_naive = time_naive.count() / time_fft.count();
        const double speedup_dir   = time_direct.count() / time_fft.count();

        std::cout << std::left
                  << std::setw(7) << order
                  << std::fixed << std::setprecision(3)
                  << std::setw(10) << time_naive.count()  << " ms "
                  << std::setw(10) << time_direct.count() << " ms "
                  << std::setw(10) << time_fft.count()    << " ms "
                  << std::setw(10) << speedup_naive       << "x  "
                  << std::setw(10) << speedup_dir         << "x  "
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << max_err_dir
                  << std::setw(15) << max_err_nai
                  << "\n" << std::flush;
    }

    // ========================================================================
    // Table 2: Isolated Azimuthal Projection (Eq. 37) Scaling & Crossover
    //   Direct O(N_theta * N_phi * M) vs FFT O(N_theta * N_phi * log N_phi)
    // ========================================================================
    std::cout << "\n========================================================================================================\n";
    std::cout << "      Optimization #4: Isolated Azimuthal Projection (Eq. 37) - Direct vs FFT Scaling\n";
    std::cout << "========================================================================================================\n";
    std::cout << "Tests both exact nPhi = 3N-2 (Bluestein) and padded nPhi = next_pow2(3N-2) (Radix-2)\n\n";

    std::cout << std::left
              << std::setw(7)  << "Order"
              << std::setw(12) << "nPhi(exact)"
              << std::setw(12) << "nPhi(pow2)"
              << std::setw(14) << "Direct"
              << std::setw(14) << "FFT(Bluest.)"
              << std::setw(14) << "FFT(Radix2)"
              << std::setw(13) << "Speedup(R2)"
              << std::setw(14) << "Max Error"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------------------------------\n" << std::flush;

    std::vector<int> az_orders = {4, 8, 12, 16, 24, 32, 48, 64};

    for (int order : az_orders) {
        const int full_order = 2 * order - 1;
        const int tt = full_order + order - 1; // 3*order - 2
        const int nTheta = (tt + 1) / 2;
        const int nPhi_exact = tt;
        const int nPhi_pow2  = sylph::AzimuthalFFT::next_pow2(tt);

        sylph::gauss_legendre gl(nTheta);
        sylph::Spherical_grid grid_exact(nTheta, nPhi_exact, gl.nodes());
        sylph::Spherical_grid grid_pow2(nTheta, nPhi_pow2, gl.nodes());

        for (int i = 0; i < nTheta; ++i) {
            for (int j = 0; j < nPhi_exact; ++j) {
                grid_exact(i, j) = dist(gen);
            }
            for (int j = 0; j < nPhi_pow2; ++j) {
                grid_pow2(i, j) = dist(gen);
            }
        }

        std::vector<double> s_direct, s_bluestein, s_radix2, s_direct_pow2;

        // Warm up plans
        sylph::AzimuthalFFT::project_rows(grid_exact, order, s_bluestein);
        sylph::AzimuthalFFT::project_rows(grid_pow2, order, s_radix2);

        const int iters = (order <= 16) ? 20 : 5;

        auto td0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) {
            direct_azimuthal_project(grid_exact, order, s_direct);
        }
        auto td1 = std::chrono::high_resolution_clock::now();
        double ms_direct = std::chrono::duration<double, std::milli>(td1 - td0).count() / iters;

        auto tb0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) {
            sylph::AzimuthalFFT::project_rows(grid_exact, order, s_bluestein);
        }
        auto tb1 = std::chrono::high_resolution_clock::now();
        double ms_bluestein = std::chrono::duration<double, std::milli>(tb1 - tb0).count() / iters;

        auto tr0 = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iters; ++it) {
            sylph::AzimuthalFFT::project_rows(grid_pow2, order, s_radix2);
        }
        auto tr1 = std::chrono::high_resolution_clock::now();
        double ms_radix2 = std::chrono::duration<double, std::milli>(tr1 - tr0).count() / iters;

        // Verify both Bluestein (on grid_exact) and Radix-2 (on grid_pow2) against direct
        direct_azimuthal_project(grid_pow2, order, s_direct_pow2);

        double max_err = 0.0;
        for (std::size_t idx = 0; idx < s_direct.size(); ++idx) {
            max_err = std::max(max_err, std::abs(s_direct[idx] - s_bluestein[idx]));
            max_err = std::max(max_err, std::abs(s_direct_pow2[idx] - s_radix2[idx]));
        }

        const double speedup_r2 = ms_direct / ms_radix2;

        std::cout << std::left
                  << std::setw(7)  << order
                  << std::setw(12) << nPhi_exact
                  << std::setw(12) << nPhi_pow2
                  << std::fixed << std::setprecision(4)
                  << std::setw(10) << ms_direct    << " ms "
                  << std::setw(10) << ms_bluestein << " ms "
                  << std::setw(10) << ms_radix2    << " ms "
                  << std::setprecision(2)
                  << std::setw(10) << speedup_r2   << "x  "
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << max_err
                  << "\n" << std::flush;
    }

    std::cout << "========================================================================================================\n\n" << std::flush;
}

} // namespace examples
