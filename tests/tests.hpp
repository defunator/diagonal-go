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

double TestColvilleFunction(const std::array<double, 4>& x) {
    return 
            + 100. * std::pow(x[0] * x[0] - x[1], 2)
            + std::pow(x[2] - 1, 2)
            + std::pow(x[0] - 1, 2)
            + 90 * std::pow(x[2] * x[2] - x[3], 2)
            + 10.1 * std::pow(x[1] - 1, 2)
            + 10.1 * std::pow(x[3] - 1, 2)
            + 19.8 * (x[3] - 1) * (x[1] - 1);
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

double TestHimmelBlauFunction(const std::array<double, 2>& x) {
    return std::pow(x[0] * x[0] + x[1] - 11, 2.) + std::pow(x[0] + x[1] * x[1] - 7, 2.);
}

template <std::size_t N>
double TestDeb3Function(const std::array<double, N>& x) {
    double res = 0;
    for (std::size_t i = 0; i != N; ++i) {
        res += std::pow(std::sin(5. * M_PI * (std::pow(x[i], 0.75) - 0.05)), 6.);
    }
    return -res / N;
}

double TestTridFunction(const std::array<double, 6>& x) {
    double res = 0;
    for (std::size_t i = 0; i != 6; ++i) {
        res += std::pow(x[i] - 1, 2.);
        if (i != 0) {
            res -= x[i] * x[i - 1];
        }
    }
    return res;
}

double TestLangerman5Function(const std::array<double, 5>& x) {
    int n = 5;
    std::vector<std::vector<double>> a = {
        {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
        {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
        {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
        {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
        {8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567}
    };

    std::vector<double> c = { 0.806, 0.517, 0.100, 0.908, 0.965 };

    double y = 0.0;
    for (int i = 0; i < 5; i++)
    { 
        double s = 0.0;
        for (int j = 0; j < n; j++)
            s += std::pow(x[j] - a[i][j], 2.);
        y += c[i] * std::exp((-1.0 / M_PI) * s) * std::cos(M_PI*s);
    }
    return -y;
}

double TestPavianiFunction(const std::array<double, 10>& x) {
    double res = 0;
    double prod = 1.;
    for (auto xx : x) {
        res += std::pow(std::log(xx - 2.), 2.);
        res += std::pow(std::log(10. - xx), 2.);
        prod *= xx;
    }
    return res - std::pow(prod, 0.2);
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

Test<4> GetTestColville() {
    double fMin = 0;
    std::array<double, 4> xMin;
    std::array<double, 4> leftBound;
    std::array<double, 4> rightBound;
    xMin.fill(1.);
    leftBound.fill(-2.);
    rightBound.fill(2.);
    assert(std::abs(TestColvilleFunction(xMin) - fMin) < 0.001);
    return Test<4>("Colville", fMin, xMin, leftBound, rightBound, TestColvilleFunction, 0.01, 1.1, 10.);
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

Test<2> GetTestHimmelBlau() {
    double fMin = 0;
    std::array<double, 2> xMin{3, 2};
    std::array<double, 2> leftBound{-6, -6};
    std::array<double, 2> rightBound{6, 6};
    assert(std::abs(TestHimmelBlauFunction(xMin) - fMin) < 0.001);
    return Test<2>("HimmelBlau", fMin, xMin, leftBound, rightBound, TestHimmelBlauFunction, 0.01, 1.2, 10.);
}

template <std::size_t N>
Test<N> GetTestDeb3() {
    double fMin = -1.;
    std::array<double, N> xMin;
    std::array<double, N> leftBound;
    std::array<double, N> rightBound;
    leftBound.fill(0);
    rightBound.fill(1);
    // assert(std::abs(TestHimmelBlauFunction(xMin) - fMin) < 0.001);
    return Test<N>("Deb 3", fMin, xMin, leftBound, rightBound, TestDeb3Function<N>, 0.01, 1.2, 10.);
}

Test<6> GetTestTrid() {
    double fMin = -50.;
    std::array<double, 6> xMin{6, 10, 12, 12, 10, 6};
    std::array<double, 6> leftBound;
    std::array<double, 6> rightBound;
    leftBound.fill(-20);
    rightBound.fill(20);
    assert(std::abs(TestTridFunction(xMin) - fMin) < 0.001);
    return Test<6>("Trid", fMin, xMin, leftBound, rightBound, TestTridFunction, 0.01, 1.3, 10.);
}

Test<5> GetTestLangerman5() {
    double fMin = -0.965;
    std::array<double, 5> xMin{8.074, 8.777, 3.467, 1.86301, 6.708};
    std::array<double, 5> leftBound;
    std::array<double, 5> rightBound;
    leftBound.fill(0);
    rightBound.fill(10);
    assert(std::abs(TestLangerman5Function(xMin) - fMin) < 0.001);
    return Test<5>("Langerman5", fMin, xMin, leftBound, rightBound, TestLangerman5Function, 0.01, 1.2, 10.);
}

Test<10> GetTestPaviani() {
    double fMin = -45.778;
    std::array<double, 10> xMin;
    std::array<double, 10> leftBound;
    std::array<double, 10> rightBound;
    xMin.fill(9.350266);
    leftBound.fill(2.1);
    rightBound.fill(9.9);
    assert(std::abs(TestPavianiFunction(xMin) - fMin) < 0.1);
    return Test<10>("Paviani", fMin, xMin, leftBound, rightBound, TestPavianiFunction, 0.01, 1.2, 10.);
}

}