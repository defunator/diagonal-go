#pragma once
#include "point_tree.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <utility>


namespace NParallel {

template<std::size_t N>
class Plane;

template<std::size_t N>
class Optimizer;


/**
 * @tparam N dimensionality of function
 */
template<std::size_t N>
class Fragment {
protected:
    std::array<std::shared_ptr<OneDimPointTree>, N> leftPointTreeRefs;
    std::array<std::shared_ptr<OneDimPointTree>, N> rightPointTreeRefs;
    double fLeft;
    double fRight;

private:
    // if direction of the diagonal is changed
    std::array<bool, N> isOddCentral;
    // std::array<std::string, N> code;
    // to divide by max len
    std::priority_queue<std::pair<double, std::size_t>> divideDim;
    std::array<std::size_t, N> leftZeroDivsCount;
    std::array<std::size_t, N> rightZeroDivsCount;

    void updDivideDim();

    void traversePointTree(char planeId);

protected:
    double diff;
    double R;
    std::size_t prevDivDim;
    double prevDiff;

public:
    template<std::size_t M>
    friend class Plane;

    template<std::size_t M>
    friend class Optimizer;

    Fragment(const std::array<double, N>& left_, const std::array<double, N>& right_,
            std::array<std::shared_ptr<OneDimPointTree>, N>&& leftPointTreeRefs_,
            std::array<std::shared_ptr<OneDimPointTree>, N>&& rightPointTreeRefs_,
            double fLeft_, double fRight_);
    Fragment(const Fragment& other);
    Fragment(Fragment&& other);

    void UpdR(double mu, double fLeft, double fRight);

    void Divide(char planeId);

    bool IsOddDim(std::size_t dim) const;

    // void TransformToPointCode(bool leftBorder = true);
    // void InvTransformToPointCode(bool leftBorder = true);
};

} // end namespace NParallel


// template<std::size_t N>
// void NParallel::Fragment<N>::TransformToPointCode(bool leftBorder) {
//     for (std::size_t i = 0; i != N; ++i) {
//         if (
//             (leftBorder && isOddCentral[i]) || (!leftBorder && !isOddCentral[i])
//         ) {
//             std::size_t j = code[i].size() - 1;
//             for (; j != 0 && code[i][j] == '2'; --j) {
//                 code[i][j] = '0';
//             }
//             code[i][j] = '0' + (code[i][j] - '0' + 1) % 3;
//         }
//     }
// }


// template<std::size_t N>
// void NParallel::Fragment<N>::InvTransformToPointCode(bool leftBorder) {
//     for (std::size_t i = 0; i != N; ++i) {
//         if (
//             (leftBorder && isOddCentral[i]) || (!leftBorder && !isOddCentral[i])
//         ) {
//             std::size_t j = code[i].size() - 1;
//             for (; j != 0 && code[i][j] == '0'; --j) {
//                 code[i][j] = '2';
//             }
//             code[i][j] = '0' + (code[i][j] - '0' + 2) % 3;
//         }
//     }
// }



template<std::size_t N>
NParallel::Fragment<N>::Fragment(
    const std::array<double, N>& left, const std::array<double, N>& right,
    std::array<std::shared_ptr<NParallel::OneDimPointTree>, N>&& leftPointTreeRefs_,
    std::array<std::shared_ptr<NParallel::OneDimPointTree>, N>&& rightPointTreeRefs_,
    double fLeft_, double fRight_)
    : leftPointTreeRefs(std::move(leftPointTreeRefs_))
    , rightPointTreeRefs(std::move(rightPointTreeRefs_))
    , fLeft(fLeft_)
    , fRight(fRight_) {
    isOddCentral.fill(false);
    // code.fill("0");
    leftZeroDivsCount.fill(0);
    rightZeroDivsCount.fill(0);
    for (std::size_t i = 0; i != N; ++i) {
        assert(left[i] < right[i]);
        divideDim.emplace(right[i] - left[i], i);
        diff += (right[i] - left[i]) * (right[i] - left[i]);
    }
    diff = std::sqrt(diff);
}


template<std::size_t N>
NParallel::Fragment<N>::Fragment(const NParallel::Fragment<N> &other)
    : leftPointTreeRefs(other.leftPointTreeRefs)
    , rightPointTreeRefs(other.rightPointTreeRefs)
    , fLeft(other.fLeft)
    , fRight(other.fRight)
    , isOddCentral(other.isOddCentral)
    // , code(other.code)
    , divideDim(other.divideDim)
    , leftZeroDivsCount(other.leftZeroDivsCount)
    , rightZeroDivsCount(other.rightZeroDivsCount)
    , diff(other.diff) { }


template<std::size_t N>
NParallel::Fragment<N>::Fragment(NParallel::Fragment<N> &&other)
    : leftPointTreeRefs(std::move(other.leftPointTreeRefs))
    , rightPointTreeRefs(std::move(other.rightPointTreeRefs))
    , fLeft(other.fLeft)
    , fRight(other.fRight)
    , isOddCentral(std::move(other.isOddCentral))
    // , code(std::move(other.code))
    , divideDim(std::move(other.divideDim))
    , leftZeroDivsCount(std::move(other.leftZeroDivsCount))
    , rightZeroDivsCount(std::move(other.rightZeroDivsCount))
    , diff(other.diff) { }


template<std::size_t N>
void NParallel::Fragment<N>::updDivideDim() {
    std::pair<double, std::size_t> dim = divideDim.top();
    divideDim.pop();
    dim.first /= 3.;
    divideDim.push(dim);
}


template<std::size_t N>
void NParallel::Fragment<N>::UpdR(double mu, double fLeft, double fRight) {
    R = mu * diff + (fLeft - fRight) * (fLeft - fRight) / (mu * diff) - 2 * (fLeft + fRight);
}

template<std::size_t N>
void NParallel::Fragment<N>::traversePointTree(char planeId) {
    std::shared_ptr<NParallel::OneDimPointTree> from;
    std::size_t zero_count;
    if (isOddCentral[prevDivDim]) {
        from = rightPointTreeRefs[prevDivDim];
        zero_count = rightZeroDivsCount[prevDivDim];

        if (planeId == '0' || planeId == '1') {
            if (planeId == '1') {
                from->Traverse('1', zero_count, leftPointTreeRefs[prevDivDim]);
            }
            leftZeroDivsCount[prevDivDim] = 0;
        } else {
            ++leftZeroDivsCount[prevDivDim];
        }
        if (planeId == '1' || planeId == '2') {
            if (planeId == '1') {
                from->Traverse('2', zero_count, rightPointTreeRefs[prevDivDim]);
            }
            rightZeroDivsCount[prevDivDim] = 0;
        } else {
            ++rightZeroDivsCount[prevDivDim];
        }
    } else {
        from = leftPointTreeRefs[prevDivDim];
        zero_count = leftZeroDivsCount[prevDivDim];

        if (planeId == '0' || planeId == '1') {
            if (planeId == '1') {
                from->Traverse('1', zero_count, rightPointTreeRefs[prevDivDim]);
            }
            rightZeroDivsCount[prevDivDim] = 0;
        } else {
            ++rightZeroDivsCount[prevDivDim];
        }
        if (planeId == '1' || planeId == '2') {
            if (planeId == '1') {
                from->Traverse('2', zero_count, leftPointTreeRefs[prevDivDim]);
            }
            leftZeroDivsCount[prevDivDim] = 0;
        } else {
            ++leftZeroDivsCount[prevDivDim];
        }
    }
}


template<std::size_t N>
void NParallel::Fragment<N>::Divide(char planeId) {
    auto[curDiff, divDim] = divideDim.top();
    prevDivDim = divDim;
    prevDiff = curDiff;

    if (planeId != '1' && isOddCentral[divDim]) {
        planeId = planeId == '0' ? '2' : '0'; 
    }
    traversePointTree(planeId);
    if (planeId == '1') {
        isOddCentral[divDim] = !isOddCentral[divDim];
    }
    // code[divDim].push_back(planeId);

    diff = diff * diff - curDiff * curDiff;
    curDiff /= 3.;
    diff += curDiff * curDiff;
    diff = std::sqrt(diff);
    updDivideDim();
}


template<std::size_t N>
bool NParallel::Fragment<N>::IsOddDim(std::size_t dim) const {
    return isOddCentral[dim];
}


