// #include "optimizer.h"
#include "plane.h"

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>


template <std::size_t N>
NSequential::Optimizer<N>::Optimizer() { }


template <std::size_t N>
double NSequential::Optimizer<N>::optimize(
    const std::array<double, N> &left,
    const std::array<double, N> &right,
    const std::function<double(const std::array<double, N>&)> &f,
    std::array<double, N>& optimum,
    double eps,
    double r,
    double C) {
    NSequential::Plane<N> plane(left, right, f, eps, r, C);
    double bigDiff2 = 0;
    for (std::size_t i = 0; i != N; ++i) {
        bigDiff2 += (right[i] - left[i]) * (right[i] - left[i]);
    }
    bigDiff2 *= eps * eps;
    while (plane.GetBestFragmentDiff2() > bigDiff2) {
        std::size_t bestFragmentId = plane.GetBestFragmentId();
        std::cout << "f(X*) = " << plane.GetBestPoint(optimum) << std::endl;
        for (auto el : optimum) {
            std::cout << el << ' ';
        }
        std::cout << std::endl;
        plane.DivideFragment(bestFragmentId);
        std::cout << plane.GetBestFragmentDiff2() << ' ' << bigDiff2 << std::endl;
        std::cout << std::endl;
    }
    return plane.GetBestPoint(optimum);
}