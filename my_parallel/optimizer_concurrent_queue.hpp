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
    double enqueueRatio = 0.9;
    double precision = 0.01;
    double updateLambdaDiff = 0.001;

    std::array<double, N> left;
    std::array<double, N> right;
    test.GetBounds(left, right);
    NParallel::Plane<N> plane(left, right, test.GetFunction(), r, C);
    double bigDiff = 0;
    for (std::size_t i = 0; i != N; ++i) {
        bigDiff += (right[i] - left[i]) * (right[i] - left[i]);
    }
    bigDiff = std::sqrt(bigDiff);
    bigDiff *= eps;
    // std::atomic<std::size_t> failInRowCnt = 0;

    omp_set_num_threads(1);
#pragma omp parallel shared(plane)
    {
    #pragma omp single
        {
            while (42) {
                if (plane.NFragments() > 50 && enqueueRatio * plane.NFragments() > plane.NQueued()) {
                    continue;
                }
                if (plane.PrevUpdMuDiff() > updateLambdaDiff) {
                    std::cout << plane.GetBestPoint(optimum) << '\n';
                    std::cout << " = ";
                    for (auto el : optimum) {
                        std::cout << el << ' ';
                    }
                    std::cout << std::endl;
                #pragma omp taskwait
                    plane.RecalcAllFragments();
                    // failInRowCnt = 0;
                }
                // if (failInRowCnt > 2.5 * (1 - enqueueRatio) * plane.NFragments()) {
                // #pragma omp taskwait
                //     break;
                // }
                if (test.GetAbsoluteValueDiff(plane.GetBestPoint()) < precision) {
                #pragma omp taskwait
                    break;
                }
                #pragma omp task
                {
                    typename NParallel::Plane<N>::NodePtr bestFragment;
                    if (plane.GetBestFragment(bestFragment)) {
                        // if (bestFragment->value->diff > bigDiff) {
                            // failInRowCnt = 0;
                            plane.DivideFragment(bestFragment);
                        // } else {
                        //     ++failInRowCnt;
                        //     plane.PutFragmentBack(bestFragment);
                        // }
                    }
                } // omp task
            }
        } // omp single
        #pragma omp taskwait
    } // omp parallel
    std::cout << "NFragments = " << plane.NFragments() << std::endl;
    std::cout << "FCount = " << plane.FCount() << std::endl;
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
