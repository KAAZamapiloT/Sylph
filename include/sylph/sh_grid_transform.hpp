#pragma once

#include "sylph/sh_coefficients.hpp"
#include "sylph/spherical_harmonics.hpp"
#include "sylph/spherical_grid.hpp"
#include "sylph/legendre.hpp"
#include "sylph/gauss_legendre.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include <types/constants.hpp>
namespace sylph{

class SHGridTransform{
public:


// equation 19
Spherical_grid  forward(
   const  SHCoefficients& coefficients,
    Spherical_grid grid
) const
{
    for (int i = 0; i < grid.theta_size(); ++i) {
        for (int j = 0; j < grid.phi_size(); ++j) {

            grid(i, j) =
                sh_.reconstruct(
                    coefficients,
                    grid.theta(i),
                    grid.phi(j)
                );
        }
    }

    return grid;
}


//eq 20 -25
Spherical_grid forward_separable(
    const SHCoefficients& coefficients,
    Spherical_grid grid
) const {

        const int order=coefficients.order();

        const int mCount=2*order-1;
          std::vector<double> s(
        static_cast<std::size_t>(grid.theta_size()) *
        static_cast<std::size_t>(mCount),
        0.0
    );

    auto mIndex=[order](int m){
        return m+order-1;
    };

     for (int i = 0; i < grid.theta_size(); ++i) {

        const double theta = grid.theta(i);

        for (int m = -order + 1;
             m <= order - 1;
             ++m)
        {
            double sum = 0.0;

            for (int l = std::abs(m);
                 l < order;
                 ++l)
            {
                sum +=
                    coefficients(l, m) *
                    C(l, m, theta);
            }

            s[
                static_cast<std::size_t>(i) * mCount
                + mIndex(m)
            ] = sum;
        }
    }


    for (int i = 0; i < grid.theta_size(); ++i) {

        for (int j = 0; j < grid.phi_size(); ++j) {

            const double phi = grid.phi(j);

            double value = 0.0;

            for (int m = -order + 1;
                 m <= order - 1;
                 ++m)
            {
                value +=
                    R(m, phi) *
                    s[
                        static_cast<std::size_t>(i) * mCount
                        + mIndex(m)
                    ];
            }

            grid(i, j) = value;
        }
    }

    return grid;

}

    SHCoefficients inverse(
        const Spherical_grid& grid,
        const std::vector<double>& weights,
        int input_order,
        int output_order
    ) const{
return inverse_impl(
        grid,
        weights,
        input_order,
        output_order,
        true
    );
}

Spherical_grid product_grid(
    std::span<const SHCoefficients* const> inputs,
    const Spherical_grid& grid
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product_grid requires at least two inputs"
        );
    }

    for (const SHCoefficients* input : inputs)
    {
        if (input == nullptr)
        {
            throw std::invalid_argument(
                "product_grid received null input"
            );
        }

        if (input->order() <= 0)
        {
            throw std::invalid_argument(
                "all input orders must be positive"
            );
        }
    }

    Spherical_grid result =
        forward_separable(*inputs[0], grid);


    for (std::size_t i = 1; i < inputs.size(); ++i)
    {
        const Spherical_grid next =
            forward_separable(*inputs[i], grid);

        result.multiply_inplace(next);
    }

    return result;
}

SHCoefficients product(
    std::span<const SHCoefficients* const> inputs,
    int output_order
) const{
    
        if (inputs.size() < 2)
        {
            throw std::invalid_argument(
                "product requires at least "
                "two inputs"
            );
        }



        int full_order = 1;

        for (const SHCoefficients* input :
             inputs)
        {
            if (input == nullptr)
            {
                throw std::invalid_argument(
                    "product received null input"
                );
            }

            if (input->order() <= 0)
            {
                throw std::invalid_argument(
                    "all input orders must be positive"
                );
            }

            full_order +=
                input->order() - 1;
        }

    
        if (output_order <= 0 ||
            output_order > full_order)
        {
            throw std::invalid_argument(
                "output_order must satisfy "
                "1 <= output_order <= full product order"
            );
        }
        const int tt=full_order + output_order-1;
        const int nTheta = (tt + 1) / 2;
        const int nPhi   = tt;

    

        gauss_legendre gl(
            nTheta
        );

        Spherical_grid grid(
            nTheta,
            nPhi,
            gl.nodes()
        );

        const Spherical_grid productGrid =
            product_grid(
                inputs,
                grid
            );

        return inverse(
            productGrid,
            gl.weights(),
            full_order,
            output_order
        );
}

SHCoefficients product(
    std::span<const SHCoefficients* const> inputs,
    int output_order,
    int t
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product requires at least two inputs"
        );
    }

    if (t <= 0)
    {
        throw std::invalid_argument(
            "t must be positive"
        );
    }



    int full_order = 1;

    for (const SHCoefficients* input : inputs)
    {
        if (input == nullptr)
        {
            throw std::invalid_argument(
                "product received null input"
            );
        }

        if (input->order() <= 0)
        {
            throw std::invalid_argument(
                "all input orders must be positive"
            );
        }

        full_order += input->order() - 1;
    }

    if (output_order <= 0 ||
        output_order > full_order)
    {
        throw std::invalid_argument(
            "output_order must satisfy "
            "1 <= output_order <= full product order"
        );
    }



    const int nTheta =
        (t + 1) / 2;

    const int nPhi =
        t;



    gauss_legendre gl(nTheta);

    Spherical_grid grid(
        nTheta,
        nPhi,
        gl.nodes()
    );

   

    const Spherical_grid productGrid =
        product_grid(
            inputs,
            grid
        );


    return inverse_impl(
        productGrid,
        gl.weights(),
        full_order,
        output_order,
        false
    );
}

double product_integral(
    std::span<const SHCoefficients* const> inputs
) const
{
    if (inputs.size() < 2)
    {
        throw std::invalid_argument(
            "product_integral requires at least "
            "two inputs"
        );
    }

    const SHCoefficients result =
        product(
            inputs,
            1
        );

    return
        2.0 *
        std::sqrt(pi) *
        result(0, 0);
}

double product_integral(
    const SHCoefficients& a,
    const SHCoefficients& b
) const
{
    const std::array<
        const SHCoefficients*,
        2
    > inputs = {
        &a,
        &b
    };

    return product_integral(
        std::span<
            const SHCoefficients* const
        >(inputs)
    );
}

 SHCoefficients product(
        const SHCoefficients& a,
        const SHCoefficients& b,
        int output_order
    ) const
    {
        const std::array<
            const SHCoefficients*,
            2
        > inputs = {
            &a,
            &b
        };

        return product(
            std::span<
                const SHCoefficients* const
            >(inputs),
            output_order
        );
    }
SHCoefficients product(
    const SHCoefficients& a,
    const SHCoefficients& b,
    int output_order,
    int t
) const
{
    const std::array<const SHCoefficients*, 2> inputs = {
        &a, &b
    };

    return product(
        std::span<const SHCoefficients* const>(inputs),
        output_order,
        t
    );
}

private:
SphericalHarmonics sh_;
Legendre  L_;

 double C(
        int l,
        int m,
        double theta
    ) 
  const  {
        const int am = std::abs(m);

        return L_.normalization(l, am) *
               L_.associated(
                   l,am,std::cos(theta)
               );
    }

    double R(
        int m,
        double phi
    ) 
   const  {
        if (m == 0) {
            return 1.0;
        }

        if (m > 0) {
            return std::sqrt(2.0) *
                   std::cos(m * phi);
        }

        return std::sqrt(2.0) *
               std::sin(std::abs(m) * phi);
    }
SHCoefficients inverse_impl(
    const Spherical_grid& grid,
    const std::vector<double>& weights,
    int input_order,
    int output_order,
    bool validate_resolution
) const
{
    if (input_order <= 0)
    {
        throw std::invalid_argument(
            "input_order must be positive"
        );
    }

    if (output_order <= 0 ||
        output_order > input_order)
    {
        throw std::invalid_argument(
            "output_order must satisfy "
            "1 <= output_order <= input_order"
        );
    }

    if (weights.size() !=
        static_cast<std::size_t>(
            grid.theta_size()
        ))
    {
        throw std::invalid_argument(
            "weights size must equal "
            "grid theta size"
        );
    }

    // ------------------------------------------------------------
    // Exact-resolution validation is optional here.
    //
    // Public inverse() passes true.
    // Approximate product(..., t) passes false.
    // ------------------------------------------------------------

    if (validate_resolution)
    {
        const int min_theta =
            (input_order + output_order) / 2;

        const int min_phi =
            input_order + output_order - 1;

        if (grid.theta_size() < min_theta)
        {
            throw std::invalid_argument(
                "theta resolution is insufficient "
                "for inverse transform"
            );
        }

        if (grid.phi_size() < min_phi)
        {
            throw std::invalid_argument(
                "phi resolution is insufficient "
                "for inverse transform"
            );
        }
    }

    const int mCount =
        2 * input_order - 1;

    auto mIndex =
        [input_order](int m)
        {
            return m + input_order - 1;
        };

    std::vector<double> s(
        static_cast<std::size_t>(
            grid.theta_size()
        ) *
        static_cast<std::size_t>(
            mCount
        ),
        0.0
    );

    const double phiOrder =
        2.0 * pi /
        static_cast<double>(
            grid.phi_size()
        );

    // ------------------------------------------------------------
    // Eq. 37: azimuthal projection
    // ------------------------------------------------------------

    for (int i = 0;
         i < grid.theta_size();
         ++i)
    {
        for (int m = -input_order + 1;
             m <= input_order - 1;
             ++m)
        {
            double sum = 0.0;

            for (int j = 0;
                 j < grid.phi_size();
                 ++j)
            {
                sum +=
                    grid(i, j) *
                    R(
                        m,
                        grid.phi(j)
                    );
            }

            s[
                static_cast<std::size_t>(i) *
                static_cast<std::size_t>(mCount) +
                static_cast<std::size_t>(
                    mIndex(m)
                )
            ] =
                phiOrder * sum;
        }
    }

    // ------------------------------------------------------------
    // Eq. 30: polar projection
    // ------------------------------------------------------------

    SHCoefficients result(
        output_order
    );

    for (int l = 0;
         l < output_order;
         ++l)
    {
        for (int m = -l;
             m <= l;
             ++m)
        {
            double coefficient = 0.0;

            for (int i = 0;
                 i < grid.theta_size();
                 ++i)
            {
                coefficient +=
                    weights[i] *
                    C(
                        l,
                        m,
                        grid.theta(i)
                    ) *
                    s[
                        static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(mCount) +
                        static_cast<std::size_t>(
                            mIndex(m)
                        )
                    ];
            }

            result(l, m) =
                coefficient;
        }
    }

    return result;
}
};

}
