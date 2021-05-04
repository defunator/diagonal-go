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


namespace NParallel {

template<std::size_t N>
class Optimizer {
private:

public:
    Optimizer();

    double optimize(
        const NTest::Test<N>& test,
        std::array<double, N>& optimum
    );

    double getScore(double fSearch, double fOptim) {
        return std::abs(fSearch - fOptim) / std::abs(fOptim);
    }
};

} // end namespace NParallel


template <std::size_t N>
NParallel::Optimizer<N>::Optimizer() { }


template <std::size_t N>
double NParallel::Optimizer<N>::optimize(
    const NTest::Test<N>& test,
    std::array<double, N>& optimum
) {
    std::array<double, N> left;
    std::array<double, N> right;
    test.GetBounds(left, right);

    double fGlobalOptimum = test.GetFunction()(left);
    optimum = left;
    std::array<double, N> lleft = left;
    std::array<double, N> rright = right;
#pragma omp parallel num_threads(8)
    {
    #pragma omp single
        {
            std::size_t nThreads = omp_get_num_threads();
            std::cout << nThreads << std::endl;
            for (std::size_t i = 1; i != nThreads + 1; ++i) {
                    lleft[0] = (left[0] * i + right[0] * (nThreads - i)) / nThreads;
                    rright[0] = (left[0] * (i - 1) + right[0] * (nThreads - i + 1)) / nThreads;
                #pragma omp task firstprivate(lleft, rright)
                    {
                        NSequential::Optimizer<N> op;
                        std::array<double, N> res;
                        double fOptimum = op.optimize(test, res);
                    #pragma omp critical
                        // std::cout << lleft[0] << ' ' << rright[0] << std::endl;
                        // std::cout << fOptimum << std::endl;
                        if (fOptimum < fGlobalOptimum) {
                            fGlobalOptimum = fOptimum;
                            optimum = res;
                        }
                    }
            }
        }
    }
    return fGlobalOptimum;
}
