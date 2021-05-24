#pragma once
#include "plane.hpp"
#include "tests.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <omp.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include <chrono>


namespace NParallel {

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

} // end namespace NParallel


template <std::size_t N>
NParallel::Optimizer<N>::Optimizer() { }

using namespace std;

template <std::size_t N>
double NParallel::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum,
    double eps,
    double r,
    double C) {
    // double enqueueRatio = 0.9;
    double precision = -43.7457- -45.778;
    double updateLambdaDiff = 0.05;

    std::array<double, N> left;
    std::array<double, N> right;
    test.GetBounds(left, right);
    NParallel::Plane<N> plane(left, right, test.GetFunction(), r, C, 1);
    double bigDiff = 0;
    for (std::size_t i = 0; i != N; ++i) {
        bigDiff += (right[i] - left[i]) * (right[i] - left[i]);
    }
    bigDiff = std::sqrt(bigDiff);
    bigDiff *= eps;
    // std::cout << bigDiff << '\n';
    // std::atomic<std::size_t> failInRowCnt = 0;

    std::size_t nThreads = 5;
    std::size_t prevUpdNFragments = 0;
    std::atomic<std::size_t> enqueued = 0;
    double bestValue = plane.GetBestPoint();
    // std::size_t maxIter = 20;
    // std::size_t curIter = 0;
    omp_set_num_threads(nThreads);
    
#pragma omp parallel shared(plane)
    {
        typename NParallel::Plane<N>::NodePtr bestFragment;
        while (42) {
            // if (plane.NFragments() > 50 && enqueueRatio * plane.NFragments() > plane.NQueued()) {
            //     continue;
            // }
            // while (enqueued >= 1. * nThreads);
            // std::cout << plane.NFragments() << ' ' << enqueued << '\n';
            if (plane.PrevUpdMuDiff() > updateLambdaDiff || plane.NFragments() - prevUpdNFragments > 0.1 * prevUpdNFragments) {
            // #pragma omp taskwait
            #pragma omp barrier
                // if (bestValue > plane.GetBestPoint()) {
                //     bestValue = plane.GetBestPoint();
                //     curIter = 0;
                // } else {
                //     ++curIter;
                // }
                // if (curIter == maxIter) {
                //     break;
                // }
                // enqueueRatio *= enqueueDecreaseRatio;
                // std::cout << enqueueRatio << '\n';
            #pragma omp single
                {
                    std::cout << plane.GetBestPoint(optimum);
                    std::cout << " = ";
                    for (auto el : optimum) {
                        std::cout << el << ' ';
                    }
                    std::cout << std::endl;
                    plane.RecalcAllFragments();
                    prevUpdNFragments = plane.NFragments();
                    // failInRowCnt = 0;
                }
            }
            // if (failInRowCnt > nThreads * 10) {
            // #pragma omp taskwait
            //     break;
            // }
            if (test.GetAbsoluteValueDiff(plane.GetBestPoint()) < precision) {
            #pragma omp barrier
                break;
            }
            // ++enqueued;
            if (plane.GetBestFragment(bestFragment)) {
                // std::cout << bestFragment->value->diff << '\n';
                // if (bestFragment->value->diff > bigDiff) {
                    // failInRowCnt = 0;
                    plane.DivideFragment(bestFragment);
                // } else {
                    // ++failInRowCnt;
                //     plane.PutFragmentBack(bestFragment);
                // }
            }
            // --enqueued;
        }
    } // omp parallel
    std::cout << "\n{\"NFragments\": " << plane.NFragments() << ", ";
    std::cout << "\"FCount\": " << plane.FCount() << ", ";
    // double res = plane.GetBestPoint(optimum);
    // std::cout << "f( ";
    // for (auto el : optimum) {
    //     std::cout << el << ' ';
    // }
    // std::cout << ") = " << res << std::endl;
    return plane.GetBestPoint(optimum);
}

template <std::size_t N>
double NParallel::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum) {
    return optimize(test, optimum, test.eps, test.r, test.C);
}
