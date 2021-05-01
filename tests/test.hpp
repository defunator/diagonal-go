#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>


namespace NTest {

template <std::size_t N>
class Test {
protected:
    std::string name;
    double fMin;
    std::array<double, N> xMin;
    std::array<double, N> leftBound;
    std::array<double, N> rightBound;
    std::function<double(const std::array<double, N>&)> f;

public:
    double eps = 0.02;
    double r = 1.2;
    double C = 10;

    Test(const std::string& name_, double fMin_, const std::array<double, N>& xMin_,
        const std::array<double, N>& leftBound_, const std::array<double, N>& rightBound_,
        const std::function<double(const std::array<double, N>&)>& f_,
        double eps_ = 0.02, double r_ = 1.2, double C_ = 10)
        : name(name_)
        , fMin(fMin_)
        , xMin(xMin_)
        , leftBound(leftBound_)
        , rightBound(rightBound_)
        , f(f_)
        , eps(eps_)
        , r(r_)
        , C(C_) {
    }

    Test(const Test& other)
        : name(other.name)
        , fMin(other.fMin)
        , xMin(other.xMin)
        , leftBound(other.leftBound)
        , rightBound(other.rightBound)
        , f(other.f)
        , eps(other.eps)
        , r(other.r)
        , C(other.C) {
    }

    const std::string& GetName() const {
        return name;
    }

    const std::array<double, N>& GetArgMin() const {
        return xMin;
    }

    std::function<double(const std::array<double, N>&)> GetFunction() const {
        return f;
    }

    double GetMinValue() const {
        return fMin;
    }

    void GetMinFragment(std::array<double, N>& to) const {
        to = xMin;
    }

    void GetBounds(std::array<double, N>& leftTo, std::array<double, N>& rightTo) const {
        leftTo = leftBound;
        rightTo = rightBound;
    }

    double GetRelativeValueDiff(double fValue) const {
        assert(fMin != 0);
        return std::abs(fMin - fValue) / std::abs(fMin);
    }

    double GetAbsoluteValueDiff(double fValue) const {
        return std::abs(fMin - fValue);
    }

    double GetRelativeArgDiff(const std::array<double, N>& x) const {
        double d2Diff = 0;
        double d2XMin = 0;
        for (std::size_t i = 0; i != N; ++i) {
            d2Diff = (x[i] - xMin[i]) * (x[i] - xMin[i]);
            d2XMin = xMin[i] * xMin[i];
        }
        assert(d2XMin != 0);
        return std::sqrt(d2Diff / d2XMin);
    }

    double GetAbsoluteArgDiff(const std::array<double, N>& x) const {
        double d2Diff = 0;
        for (std::size_t i = 0; i != N; ++i) {
            d2Diff += (x[i] - xMin[i]) * (x[i] - xMin[i]);
        }
        return d2Diff;
    }
}; // end Test

} // end NTest
