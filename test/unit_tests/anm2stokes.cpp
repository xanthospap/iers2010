#include "stokes_coefficients.hpp"
#include "eigen3/Eigen/Eigen"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

// @snippet example_usage

/** Tests to validate the StokesCoeffs.from_anm() function.
 */

std::uniform_real_distribution<double> unif(-1e0, 1e0);
std::default_random_engine re;

using namespace dso;

int main() {

    // Test 1: Small square case (deg == ord)
    {
        printf("Test Case 1: Small Square Matrix\n");
        const int deg = 10;
        const int ord = 10;
        StokesCoeffs s(deg);
        
        // Populate C and S with values
        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = 1e0;
                if (m!=0) s.S(l, m) = 2e0;
            }
        }

        auto anm = s.to_anm();

        // Create a new StokesCoeffs object and reconstruct C and S
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm, ord);

        // Verify values
        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= row; col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (col != 0) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == ord);
    }
    
    // Test 2: Small square case (deg == ord)
    {
        printf("Test Case 2: Small Square Matrix\n");
        const int deg = 10;
        const int ord = 10;
        StokesCoeffs s(deg);
        
        // Populate C and S with values
        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = 1e0;
                if (m!=0) s.S(l, m) = 2e0;
            }
        }

        auto anm = s.to_anm();

        // Create a new StokesCoeffs object and reconstruct C and S
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm);

        // Verify values
        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= row; col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (0 != col) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == ord);
    }

    // Test 3: Randomized values, square case
    {
        printf("Test Case 3: Randomized Square Matrix\n");
        const int deg = 20;
        const int ord = 20;
        StokesCoeffs s(deg);

        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = unif(re);
                if (m!=0) s.S(l, m) = unif(re);
            }
        }

        auto anm = s.to_anm();
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm, ord);

        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= row; col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (0 != col) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == ord);
    }

    // Test 4: Trapezoidal case (deg != ord)
    {
        printf("Test Case 4: Trapezoidal Matrix (deg != ord)\n");
        const int deg = 26;
        const int ord = 19;
        StokesCoeffs s(deg, ord);

        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = unif(re);
                if (m!=0) s.S(l, m) = unif(re);
            }
        }

        auto anm = s.to_anm();
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm, ord);

        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= std::min(row, ord); col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (0!= col) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == ord);
    }
    
    // Test 5: Trapezoidal case (deg != ord)
    {
        printf("Test Case 5: Trapezoidal Matrix (deg != ord)\n");
        const int deg = 126;
        const int ord = 125;
        StokesCoeffs s(deg, ord);

        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = unif(re);
                if (m!=0) s.S(l, m) = unif(re);
            }
        }

        auto anm = s.to_anm();
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm, ord);

        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= std::min(row, ord); col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (0 != col) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == ord);
    }
    
    // Test 6: Trapezoidal case (deg != ord)
    {
        printf("Test Case 4: Trapezoidal Matrix (deg != ord)\n");
        const int deg = 126;
        const int ord = 125;
        StokesCoeffs s(deg, ord);

        for (int m = 0; m <= ord; m++) {
            for (int l = m; l <= deg; l++) {
                s.C(l, m) = unif(re);
                if (m!=0) s.S(l, m) = unif(re);
            }
        }

        auto anm = s.to_anm();
        StokesCoeffs reconstructed_s = StokesCoeffs::from_anm(anm);

        for (int row = 0; row <= deg; row++) {
            for (int col = 0; col <= std::min(row, ord); col++) {
                assert(std::abs(reconstructed_s.C(row, col) - s.C(row, col)) < 1e-10);
                if (0 != col) {
                    assert(std::abs(reconstructed_s.S(row, col) - s.S(row, col)) < 1e-10);
                }
            }
        }
        
        assert(reconstructed_s.max_degree() == deg);
        assert(reconstructed_s.max_order() == deg);
    }
    
    return 0;
}
