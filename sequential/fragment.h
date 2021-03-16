#pragma once
#include <array>
#include <iostream>
#include <set>
#include <string>
#include <utility>


namespace NSequential {
/**
 * @tparam N dimensionality of function
 */
template<std::size_t N>
class Fragment {
public:
    std::array<std::string, N> code;
    std::string isOddCentral;
    std::set<std::pair<double, std::size_t>> divideDim;
    double diff2;
    double R;

    Fragment(const std::array<std::string, N>& code_,
            const std::set<std::pair<double, std::size_t>>& divideDim_);
    Fragment(std::array<std::string, N>&& code_,
            std::set<std::pair<double, std::size_t>>&& divideDim_);
    Fragment(const Fragment& other);
    Fragment(Fragment&& other);

    void TransformToPointCode(bool leftBorder = true);
    void InvTransformToPointCode(bool leftBorder = true);

    std::size_t getDivideDim();
    void updDivideDim();

    void updR(double C, double r, std::size_t k, double dMax,
                double lambdaMax, double fLeft, double fRight);

    void Divide(char planeId, Fragment &other);
    void InplaceDivide(char planeId);
};

}

#include "fragment.tpp"