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
    // if direction of the diagonal is changed
    std::array<bool, N> isOddCentral;
    // to divide by max len
    std::set<std::pair<double, std::size_t>> divideDim;

    void updDivideDim();

public:
    double diff;
    double R;

    Fragment(const std::array<std::string, N>& code_,
            const std::set<std::pair<double, std::size_t>>& divideDim_);
    Fragment(std::array<std::string, N>&& code_,
            std::set<std::pair<double, std::size_t>>&& divideDim_);
    Fragment(const Fragment& other);
    Fragment(Fragment&& other);

    void TransformToPointCode(bool leftBorder = true);
    void InvTransformToPointCode(bool leftBorder = true);

    std::pair<double, std::size_t> GetDivideDim();

    void UpdR(double mu, double fLeft, double fRight);

    void Divide(char planeId);

    const std::array<std::string, N>& GetCode() const;

    bool IsOddDim(std::size_t dim) const;
};

} // end namespace NSequential


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    const std::array<std::string, N>& code_,
    const std::set<std::pair<double, std::size_t>>& divideDim_)
    : code(code_)
    , divideDim(divideDim_) {
    isOddCentral.fill(false);
    diff = 0;    
    for (auto& el : divideDim) {
        diff += el.first * el.first;
    }
    diff = std::sqrt(diff);
}


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    std::array<std::string, N>&& code_,
    std::set<std::pair<double, std::size_t>>&& divideDim_)
    : code(std::move(code_))
    , divideDim(std::move(divideDim_)) {
    isOddCentral.fill(false);
    diff = 0;    
    for (auto& el : divideDim) {
        diff += el.first * el.first;
    }
    diff = std::sqrt(diff);
}


template<std::size_t N>
NSequential::Fragment<N>::Fragment(const NSequential::Fragment<N> &other)
    : code(other.code)
    , isOddCentral(other.isOddCentral)
    , divideDim(other.divideDim)
    , diff(other.diff) { }


template<std::size_t N>
NSequential::Fragment<N>::Fragment(NSequential::Fragment<N> &&other)
    : code(std::move(other.code))
    , isOddCentral(std::move(other.isOddCentral))
    , divideDim(std::move(other.divideDim))
    , diff(other.diff) { }


template<std::size_t N>
void NSequential::Fragment<N>::TransformToPointCode(bool leftBorder) {
    for (std::size_t i = 0; i != N; ++i) {
        // diagonal changed direction
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
std::pair<double, std::size_t> NSequential::Fragment<N>::GetDivideDim() {
    return *divideDim.rbegin();
}


template<std::size_t N>
void NSequential::Fragment<N>::updDivideDim() {
    std::pair<double, std::size_t> dim = *divideDim.rbegin();
    divideDim.erase(--divideDim.end());
    dim.first /= 3.;
    divideDim.insert(dim);
}


template<std::size_t N>
void NSequential::Fragment<N>::UpdR(double mu, double fLeft, double fRight) {
    R = mu * diff + (fLeft - fRight) * (fLeft - fRight) / (mu * diff) - 2 * (fLeft + fRight);
}


template<std::size_t N>
void NSequential::Fragment<N>::Divide(char planeId) {
    auto[curDiff, divDim] = GetDivideDim();
    code[divDim].push_back(planeId);
    if (planeId == '1') {
        isOddCentral[divDim] = !isOddCentral[divDim];
    }
    diff = diff * diff - curDiff * curDiff;
    curDiff /= 3.;
    diff += curDiff * curDiff;
    diff = std::sqrt(diff);
    updDivideDim();
}


template<std::size_t N>
const std::array<std::string, N>& NSequential::Fragment<N>::GetCode() const {
    return code;
}


template<std::size_t N>
bool NSequential::Fragment<N>::IsOddDim(std::size_t dim) const {
    return isOddCentral[dim];
}


