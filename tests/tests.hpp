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

double TestChichinadzeFunction(const std::array<double, 2>& x) {
    return x[0] * x[0] - 12. * x[0] + 11. + 10. * std::cos(0.5 * M_PI * x[0])
            + 8 * std::sin(2.5 * M_PI * x[0]) - std::sqrt(0.2) * std::exp(-0.5 * std::pow(x[1] - 0.5, 2));
}

double TestDolanFunction(const std::array<double, 5>& x) {
    return std::abs((x[0] + 1.7 * x[1]) * std::sin(x[0]) - 1.5 * x[2]
                - 0.1 * x[3] * std::cos(x[3] + x[4] - x[0]) + 0.2 * x[4] * x[4]
                - x[1] - 1);
}

double TestHartmanFunction(const std::array<double, 6>& x) {
    std::vector<std::vector<double>> a = { { 10, 3, 17, 3.5, 1.7, 8 } ,
                                            { 0.05, 10, 17, 0.1, 8, 14 },
                                            { 3, 3.5, 1.7, 10, 17, 8 }, 
                                            { 17, 8, 0.05, 10, 0.1, 14 } };
    std::vector<std::vector<double>> p = { { 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886 },
                                            { 0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991 },
                                            { 0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650 },
                                            { 0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381 } };
    std::vector<double> c = { 1.0, 1.2, 3.0, 3.2 };

    double y = 0.0;
    for (int i = 0; i < 4; i++)
    {
        double e = 0.0;
        for (int j = 0; j < 6; j++)
            e += a[i][j] * std::pow(x[j] - p[i][j], 2.);
        y += c[i] * std::exp(-e);
    }
    return -y;
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


Test<2> GetTestChichinadze() {
    double fMin = -42.94438701899098;
    std::array<double, 2> xMin{6.189866586965680, 0.5};
    std::array<double, 2> leftBound;
    std::array<double, 2> rightBound;
    leftBound.fill(-30.);
    rightBound.fill(30.);
    assert(std::abs(TestChichinadzeFunction(xMin) - fMin) < 0.001);
    return Test<2>("Chichinadze", fMin, xMin, leftBound, rightBound, TestChichinadzeFunction, 0.01, 1.3, 10.);
}

Test<5> GetTestDolan() {
    double fMin = 0;
    std::array<double, 5> xMin{-74.10522498, 44.33511286, 6.21069214, 18.42772233, -16.5839403};
    std::array<double, 5> leftBound;
    std::array<double, 5> rightBound;
    leftBound.fill(-100.);
    rightBound.fill(100.);
    assert(std::abs(TestDolanFunction(xMin) - fMin) < 0.001);
    return Test<5>("Dolan", fMin, xMin, leftBound, rightBound, TestDolanFunction, 0.01, 1.1, 10.);
}

Test<6> GetTestHartman() {
    double fMin = -3.32236;
    std::array<double, 6> xMin{0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.657301};
    std::array<double, 6> leftBound;
    std::array<double, 6> rightBound;
    leftBound.fill(0.);
    rightBound.fill(1.);
    assert(std::abs(TestHartmanFunction(xMin) - fMin) < 0.001);
    return Test<6>("Hartman", fMin, xMin, leftBound, rightBound, TestHartmanFunction, 0.01, 1.1, 10.);
}

}