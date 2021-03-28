#pragma once

#include <array>
#include <cmath>
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
private:
    std::array<std::string, N> code;
    std::string isOddCentral;
    std::set<std::pair<double, std::size_t>> divideDim;

public:
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

    const std::array<std::string, N>& getCode() const;

    bool isOddDim(std::size_t dim) const;
};

} // end namespace NSequential


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    const std::array<std::string, N>& code_,
    const std::set<std::pair<double, std::size_t>>& divideDim_)
    : code(code_)
    , isOddCentral(std::string(N, false))
    , divideDim(divideDim_) {
    diff2 = 0;    
    for (auto& el : divideDim) {
        diff2 += el.first * el.first;
    }
}


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    std::array<std::string, N>&& code_,
    std::set<std::pair<double, std::size_t>>&& divideDim_)
    : code(std::move(code_))
    , isOddCentral(std::string(N, false))
    , divideDim(std::move(divideDim_)) {
    diff2 = 0;    
    for (auto& el : divideDim) {
        diff2 += el.first * el.first;
    }
}


template<std::size_t N>
NSequential::Fragment<N>::Fragment(const NSequential::Fragment<N> &other)
    : code(other.code)
    , isOddCentral(other.isOddCentral)
    , divideDim(other.divideDim)
    , diff2(other.diff2) { }


template<std::size_t N>
NSequential::Fragment<N>::Fragment(NSequential::Fragment<N> &&other)
    : code(std::move(other.code))
    , isOddCentral(std::move(other.isOddCentral))
    , divideDim(std::move(other.divideDim))
    , diff2(other.diff2) { }


template<std::size_t N>
void NSequential::Fragment<N>::TransformToPointCode(bool leftBorder) {
    for (std::size_t i = 0; i != N; ++i) {
        if (leftBorder && isOddCentral[i]) {
            std::size_t j = code[i].size() - 1;
            for (; j != 0 && code[i][j] == '2'; --j) {
                code[i][j] = '0';
            }
            code[i][j] = '0' + (code[i][j] - '0' + 1) % 3;
        } else if (!leftBorder && !isOddCentral[i]) {
            std::size_t j = code[i].size() - 1;
            for (; j != 0 && code[i][j] == '2'; --j) {
                code[i][j] = '0';
            }
            code[i][j] = '0' + (code[i][j] - '0' + 1) % 3;
        }
    }
}


template<std::size_t N>
void NSequential::Fragment<N>::InvTransformToPointCode(bool leftBorder) {
    for (std::size_t i = 0; i != N; ++i) {
        if (leftBorder && isOddCentral[i]) {
            std::size_t j = code[i].size() - 1;
            for (; j != 0 && code[i][j] == '0'; --j) {
                code[i][j] = '2';
            }
            code[i][j] = '0' + (code[i][j] - '0' + 2) % 3;
        } else if (!leftBorder && !isOddCentral[i]) {
            std::size_t j = code[i].size() - 1;
            for (; j != 0 && code[i][j] == '0'; --j) {
                code[i][j] = '2';
            }
            code[i][j] = '0' + (code[i][j] - '0' + 2) % 3;
        }
    }
}


template<std::size_t N>
std::size_t NSequential::Fragment<N>::getDivideDim() {
    return divideDim.rbegin()->second;
}


template<std::size_t N>
void NSequential::Fragment<N>::updDivideDim() {
    std::pair<double, std::size_t> dim = *divideDim.rbegin();
    divideDim.erase(--divideDim.end());
    dim.first /= 3.;
    divideDim.insert(dim);
}


template<std::size_t N>
void NSequential::Fragment<N>::updR(
    double C, double r, std::size_t k,
    double dMax, double lambdaMax, double fLeft, double fRight
) {
    double diff = std::sqrt(diff2);
    double lambda = std::abs(fRight - fLeft) / diff;
    double gamma = lambdaMax * diff / dMax;
    double mu = (r + C / k) * std::max(lambda, std::max(
        std::numeric_limits<double>::epsilon(),
        gamma
    ));
    // R = mu * diff / 2 - (fLeft + fRight) / 2;
    mu = (r + C / k) * lambdaMax;
    R = mu * diff + (fLeft - fRight) * (fLeft - fRight) / (mu * diff) - 2 * (fLeft + fRight);
}


template<std::size_t N>
void NSequential::Fragment<N>::Divide(char planeId,
                                        NSequential::Fragment<N> &other) {
    // other.code = code;
    // other.divideDim = divideDim;
    // other.diff2 = diff2;
    std::size_t divDim = getDivideDim();
    other.code[divDim].push_back(planeId);
    if (planeId == '1') {
        other.isOddCentral[divDim] = ~other.isOddCentral[divDim];
    }
    double curDiff = divideDim.rbegin()->first;
    other.diff2 -= curDiff * curDiff;
    curDiff /= 3.;
    other.diff2 += curDiff * curDiff;
    other.updDivideDim();
}


template<std::size_t N>
void NSequential::Fragment<N>::InplaceDivide(char planeId) {
    std::size_t divDim = getDivideDim();
    code[divDim].push_back(planeId);
    if (planeId == '1') {
        isOddCentral[divDim] = ~isOddCentral[divDim];
    }
    double curDiff = divideDim.rbegin()->first;
    diff2 -= curDiff * curDiff;
    curDiff /= 3.;
    diff2 += curDiff * curDiff;
    updDivideDim();
}


template<std::size_t N>
const std::array<std::string, N>& NSequential::Fragment<N>::getCode() const {
    return code;
}


template<std::size_t N>
bool NSequential::Fragment<N>::isOddDim(std::size_t dim) const {
    return isOddCentral[dim];
}

