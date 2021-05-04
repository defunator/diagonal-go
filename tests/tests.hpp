#pragma once

#include "test.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


namespace NTest {

const std::array<std::string, 3> TestNamesList{
    "Test 1",
    "Test 2",
    "Test 3"
};

double Test0Function(const std::array<double, 1>& x) {
    return x[0] * x[0];
}

double Test1Function(const std::array<double, 2>& x) {
    return 0.25 * std::pow(x[0], 4) - 0.5 * std::pow(x[0], 2) + 0.1 * x[0] + 0.5 * std::pow(x[1], 2);
}

double Test2Function(const std::array<double, 3>& x) {
    return (x[0] * x[0] - 2 * x[1] * x[1] + x[2] * x[2]) * sin(x[0]) * sin(x[1]) * sin(x[2]);
}

double Test3Function(const std::array<double, 3>& x) {
    return 1 + 100 * (x[2] - 0.25 * std::pow(x[0] + x[1], 2)) * (x[2] - 0.25 * std::pow(x[0] + x[1], 2)) + std::pow(1. - x[0], 2) + std::pow(1. - x[1], 2);
}

template <std::size_t N>
double Test4Function(const std::array<double, N>& x) {
    std::array<double, N> y;
    for (std::size_t i = 0; i != N; ++i) {
        y[i] = 1. + (x[i] - 1) / 4;
    }
    double f = 0;
    f += 10. * sin(M_PI * y[0]) * sin(M_PI * y[0]);
    f += (y[N-1] - 1) * (y[N-1] - 1);
    for (std::size_t i = 0; i != N-1; ++i) {
        f += (y[i] - 1) * (y[i] - 1) * (1 + 10 * sin(M_PI * y[i + 1]) * sin(M_PI * y[i + 1]));
    }
    return M_PI * f / 10. + 1.;
}

double Test5Function(const std::array<double, 3>& x) {
    return (std::pow(x[0], 2) - 2. * std::pow(x[1], 2) + std::pow(x[2], 2)) * std::sin(x[0]) * std::sin(x[1]) * std::sin(x[2]);
}

template <std::size_t N>
double Test6Function(const std::array<double, N>& x) {
    double f = 0;
    f += std::pow(std::sin(3 * M_PI * x[0]), 2);
    for (size_t i = 1; i != N - 1; ++i) {
        f += std::pow(x[i] - 1, 2) * (1 + std::pow(std::sin(3 * M_PI * x[i + 1]), 2));
    }
    f += std::pow(x[N-1] - 1, 2) * (1 + std::pow(std::sin(2 * M_PI * x[N-1]), 2));
    f *= 0.1;
    return f + 1;
}


Test<1> GetTest0() {
    double fMin = 0;
    std::array<double, 1> xMin{0};
    std::array<double, 1> leftBound{-10.};
    std::array<double, 1> rightBound{10.};
    return Test<1>("Test 0", fMin, xMin, leftBound, rightBound, Test0Function, 0.01, 1.1, 10.);
}

Test<2> GetTest1() {
    double fMin = -0.352386;
    std::array<double, 2> xMin{-1.04668, 0.};
    std::array<double, 2> leftBound{-10., -10.};
    std::array<double, 2> rightBound{10., 10.};
    return Test<2>("Test 1", fMin, xMin, leftBound, rightBound, Test1Function, 0.01, 1.1, 10.);
}

Test<3> GetTest2() {
    double fMin = -0.516374;
    std::array<double, 3> xMin{1., 0.556002, -1.};
    std::array<double, 3> leftBound{-0.9, -1., -1.};
    std::array<double, 3> rightBound{1., 1., 1.};
    return Test<3>("Test 2", fMin, xMin, leftBound, rightBound, Test2Function, 0.02, 1.2, 100.);
}

Test<3> GetTest3() {
    double fMin = 1.;
    std::array<double, 3> xMin{1., 1., 1.};
    std::array<double, 3> leftBound{0., 0., 0.};
    std::array<double, 3> rightBound{2., 2., 2.};
    return Test<3>("Test 3", fMin, xMin, leftBound, rightBound, Test3Function, 0.02, 1.2, 100.);
}

template <std::size_t N>
Test<N> GetTest4() {
    double fMin = 1.;
    std::array<double, N> xMin;
    std::array<double, N> leftBound;
    std::array<double, N> rightBound;
    xMin.fill(1.);
    leftBound.fill(-10.);
    rightBound.fill(10.);
    return Test<N>("Test 4", fMin, xMin, leftBound, rightBound, Test4Function<N>, 0.01, 1.1, 10.);
}

Test<3> GetTest5() {
    double fMin = -0.516374;
    std::array<double, 3> xMin{1., 0.556002, -1.};
    std::array<double, 3> leftBound;
    std::array<double, 3> rightBound;
    leftBound.fill(-1.);
    rightBound.fill(1.);
    return Test<3>("Test 5", fMin, xMin, leftBound, rightBound, Test5Function, 0.02, 1.2, 100.);
}

template <std::size_t N>
Test<N> GetTest6() {
    double fMin = 1.;
    std::array<double, N> xMin;
    std::array<double, N> leftBound;
    std::array<double, N> rightBound;
    xMin.fill(1.);
    leftBound.fill(-5.);
    rightBound.fill(5.);
    return Test<N>("Test 6", fMin, xMin, leftBound, rightBound, Test6Function<N>, 0.01, 1.1, 10.);
}



}