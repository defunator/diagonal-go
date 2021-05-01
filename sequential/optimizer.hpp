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
    while (plane.GetBestFragmentDiff() > bigDiff) {
        // std::cout << plane.GetBestFragmentDiff() << std::endl;
        // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::size_t bestFragmentId = plane.GetBestFragmentId();
        // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // std::cout << "GetPointIdByCode TIME = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
    
        // std::cout << "f(X*) = " << plane.GetBestPoint(optimum) << std::endl;
        // for (auto el : optimum) {
        //     std::cout << el << ' ';
        // }
        // std::cout << std::endl;
        // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        plane.DivideFragment(bestFragmentId);
        ++i;
        // if (i == ) break;
        // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // std::cout << "DivideFragment TIME = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
    
        // std::cout << plane.GetBestFragmentDiff() << ' ' << bigDiff << std::endl;
        // std::cout << std::endl;
        // if (test.GetRelativeValueDiff(plane.GetBestPoint(optimum)) < 0.1) {
        //     break;
        // }
        // if (i % 30 == 0) {
        //     std::cout << i << ' ' << plane.GetBestFragmentDiff() << ' ' << bigDiff << ' ' << plane.FCount() << std::endl;
        //     std::cout << plane.GetBestPoint(optimum) << '\n';
        //     for (auto x : optimum) {
        //         std::cout << x << ' ';
        //     }
        //     std::cout << '\n';
        // }
    }
    std::cout << plane.searchFragments.size() << std::endl;
    std::cout << plane.FCount() << std::endl;
    return plane.GetBestPoint(optimum);
}

template <std::size_t N>
double NSequential::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum) {
    return optimize(test, optimum, test.eps, test.r, test.C);
}
