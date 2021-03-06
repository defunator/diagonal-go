#pragma once
#include "plane.hpp"
#include "tests.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>

#include <chrono>


namespace NSequential {

/**
 * Diagonal sequential optimizer for minimization
 * @tparam N dimensionality of optimized function
 */
template <std::size_t N>
class Optimizer {
private:

public:
    Optimizer();

    /**
     * Optimization will be done on this space:
     * @f$X = \{x | \forall i = 1 \cdots n \, left[i] \leq x[i] \leq right[i]\}@f$
     * @param f Lipschitz function to optimize
     * @param left Left point of hypercube optimization space
     * @param right Right point of hypercube optimization space
     */
    double optimize(
        const NTest::Test<N>& test,
        std::array<double, N>& optimum,
        double eps,
        double r,
        double C
    );

    double optimize(
        const NTest::Test<N>& test,
        std::array<double, N>& optimum
    );
};

} // end namespace NSequential


template <std::size_t N>
NSequential::Optimizer<N>::Optimizer() { }


template <std::size_t N>
double NSequential::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum,
    double eps,
    double r,
    double C) {
    std::array<double, N> left;
    std::array<double, N> right;
    test.GetBounds(left, right);
    NSequential::Plane<N> plane(left, right, test.GetFunction(), r, C);
    double bigDiff = 0;
    for (std::size_t i = 0; i != N; ++i) {
        bigDiff += (right[i] - left[i]) * (right[i] - left[i]);
    }
    bigDiff = std::sqrt(bigDiff);
    bigDiff *= eps;
    int i = 0;
    double mean_divide = 0;
    double mean_get_best = 0;

    std::size_t bestFragmentId = plane.GetBestFragmentId();
    while (plane.GetBestFragmentDiff() > bigDiff) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        plane.DivideFragment(bestFragmentId);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        mean_divide = 1. * (1. * i * mean_divide + std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / (i + 1);

        begin = std::chrono::steady_clock::now();
        bestFragmentId = plane.GetBestFragmentId();
        end = std::chrono::steady_clock::now();
        mean_get_best = 1. * (1. * i * mean_get_best + std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / (i + 1);

        ++i;
        // if (i % 100 == 0) {
        //     std::cout << "mean_divide = " << mean_divide << std::endl;
        //     std::cout << "mean_get_best = " << mean_get_best << std::endl;
        // }
    }
    std::cout << "mean_divide = " << mean_divide << std::endl;
    std::cout << "mean_get_best = " << mean_get_best << std::endl;
    std::cout << "searchFragments.size() = " << plane.searchFragments.size() << std::endl;
    std::cout << "FCount = " << plane.FCount() << std::endl;
    return plane.GetBestPoint(optimum);
}

template <std::size_t N>
double NSequential::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum) {
    return optimize(test, optimum, test.eps, test.r, test.C);
}
