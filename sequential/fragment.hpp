#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <utility>


namespace NSequential {

template<std::size_t N>
class Plane;

/**
 * @tparam N dimensionality of function
 */
template<std::size_t N>
class Fragment {
protected:
    std::array<std::string, N> code;

private:
    // if direction of the diagonal is changed
    std::array<bool, N> isOddCentral;
    // to divide by max len
    std::priority_queue<std::pair<double, std::size_t>> divideDim;

    void updDivideDim();

public:
    double diff;
    std::size_t leftPointId;
    std::size_t rightPointId;
    double R;

    template<std::size_t M>
    friend class Plane;

    Fragment(std::array<std::string, N>&& code_,
            std::priority_queue<std::pair<double, std::size_t>>&& divideDim_,
            double diff_);
    Fragment(const Fragment& other);
    Fragment(Fragment&& other);

    void TransformToPointCode(bool leftBorder = true);
    void InvTransformToPointCode(bool leftBorder = true);

    void UpdR(double mu, double fLeft, double fRight);

    void Divide(char planeId);

    bool IsOddDim(std::size_t dim) const;
};

} // end namespace NSequential


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    std::array<std::string, N>&& code_,
    std::priority_queue<std::pair<double, std::size_t>>&& divideDim_,
    double diff_)
    : code(std::move(code_))
    , divideDim(std::move(divideDim_))
    , diff(diff_) {
    isOddCentral.fill(false);
}


template<std::size_t N>
NSequential::Fragment<N>::Fragment(const NSequential::Fragment<N> &other)
    : code(other.code)
    , isOddCentral(other.isOddCentral)
    , divideDim(other.divideDim)
    , diff(other.diff)
    , leftPointId(other.leftPointId)
    , rightPointId(other.rightPointId) { }


template<std::size_t N>
NSequential::Fragment<N>::Fragment(NSequential::Fragment<N> &&other)
    : code(std::move(other.code))
    , isOddCentral(std::move(other.isOddCentral))
    , divideDim(std::move(other.divideDim))
    , diff(other.diff)
    , rightPointId(other.rightPointId)
    , leftPointId(other.leftPointId) { }


template<std::size_t N>
void NSequential::Fragment<N>::TransformToPointCode(bool leftBorder) {
    for (std::size_t i = 0; i != N; ++i) {
        if (
            (leftBorder && isOddCentral[i]) || (!leftBorder && !isOddCentral[i])
        ) {
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
        if (
            (leftBorder && isOddCentral[i]) || (!leftBorder && !isOddCentral[i])
        ) {
            std::size_t j = code[i].size() - 1;
            for (; j != 0 && code[i][j] == '0'; --j) {
                code[i][j] = '2';
            }
            code[i][j] = '0' + (code[i][j] - '0' + 2) % 3;
        }
    }
}


template<std::size_t N>
void NSequential::Fragment<N>::updDivideDim() {
    std::pair<double, std::size_t> dim = divideDim.top();
    divideDim.pop();
    dim.first /= 3.;
    divideDim.push(dim);
}


template<std::size_t N>
void NSequential::Fragment<N>::UpdR(double mu, double fLeft, double fRight) {
    R = mu * diff + (fLeft - fRight) * (fLeft - fRight) / (mu * diff) - 2 * (fLeft + fRight);
}


template<std::size_t N>
void NSequential::Fragment<N>::Divide(char planeId) {
    auto[curDiff, divDim] = divideDim.top();
    code[divDim].push_back(planeId);
    if (planeId == '1') {
        isOddCentral[divDim] = !isOddCentral[divDim];
    } else if (isOddCentral[divDim]) {
        code[divDim].back() = planeId == '0' ? '2' : '0'; 
    }
    diff = diff * diff - curDiff * curDiff;
    curDiff /= 3.;
    diff += curDiff * curDiff;
    diff = std::sqrt(diff);
    updDivideDim();
}


template<std::size_t N>
bool NSequential::Fragment<N>::IsOddDim(std::size_t dim) const {
    return isOddCentral[dim];
}


