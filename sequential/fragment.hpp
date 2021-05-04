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


namespace NSequential {

template<std::size_t N>
class Plane;

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
    std::size_t leftPointId;
    std::size_t rightPointId;

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

    Fragment(const std::array<double, N>& left_, const std::array<double, N>& right_,
            std::array<std::shared_ptr<NSequential::OneDimPointTree>, N>&& leftPointTreeRefs_,
            std::array<std::shared_ptr<NSequential::OneDimPointTree>, N>&& rightPointTreeRefs_,
            double fLeft_, double fRight_);
    Fragment(const Fragment& other);
    Fragment(Fragment&& other);

    void TransformToPointCode(bool leftBorder = true);
    void InvTransformToPointCode(bool leftBorder = true);
    // void DimTransformToPointCode(std::size_t dim, bool leftBorder = true);
    // void DimInvTransformToPointCode(std::size_t dim, bool leftBorder = true);

    void UpdR(double mu, double fLeft, double fRight);

    void Divide(char planeId);

    bool IsOddDim(std::size_t dim) const;
};

} // end namespace NSequential


template<std::size_t N>
NSequential::Fragment<N>::Fragment(
    const std::array<double, N>& left, const std::array<double, N>& right,
    std::array<std::shared_ptr<NSequential::OneDimPointTree>, N>&& leftPointTreeRefs_,
    std::array<std::shared_ptr<NSequential::OneDimPointTree>, N>&& rightPointTreeRefs_,
    double fLeft_, double fRight_)
    : leftPointTreeRefs(std::move(leftPointTreeRefs_))
    , rightPointTreeRefs(std::move(rightPointTreeRefs_))
    , fLeft(fLeft_)
    , fRight(fRight_)
    , leftPointId(1)
    , rightPointId(2) {
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
NSequential::Fragment<N>::Fragment(const NSequential::Fragment<N> &other)
    : leftPointTreeRefs(other.leftPointTreeRefs)
    , rightPointTreeRefs(other.rightPointTreeRefs)
    , fLeft(other.fLeft)
    , fRight(other.fRight)
    , leftPointId(other.leftPointId)
    , rightPointId(other.rightPointId)
    , isOddCentral(other.isOddCentral)
    // , code(other.code)
    , divideDim(other.divideDim)
    , leftZeroDivsCount(other.leftZeroDivsCount)
    , rightZeroDivsCount(other.rightZeroDivsCount)
    , diff(other.diff) { }


template<std::size_t N>
NSequential::Fragment<N>::Fragment(NSequential::Fragment<N> &&other)
    : leftPointTreeRefs(std::move(other.leftPointTreeRefs))
    , rightPointTreeRefs(std::move(other.rightPointTreeRefs))
    , fLeft(other.fLeft)
    , fRight(other.fRight)
    , leftPointId(other.leftPointId)
    , rightPointId(other.rightPointId)
    , isOddCentral(std::move(other.isOddCentral))
    // , code(std::move(other.code))
    , divideDim(std::move(other.divideDim))
    , leftZeroDivsCount(std::move(other.leftZeroDivsCount))
    , rightZeroDivsCount(std::move(other.rightZeroDivsCount))
    , diff(other.diff) { }


// template<std::size_t N>
// void NSequential::Fragment<N>::TransformToPointCode(bool leftBorder) {
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
// void NSequential::Fragment<N>::InvTransformToPointCode(bool leftBorder) {
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


// template<std::size_t N>
// void NSequential::Fragment<N>::DimTransformToPointCode(std::size_t dim, bool leftBorder) {
//     if (
//         (leftBorder && isOddCentral[dim]) || (!leftBorder && !isOddCentral[dim])
//     ) {
//         std::size_t j = code[dim].size() - 1;
//         for (; j != 0 && code[dim][j] == '2'; --j) {
//             code[dim][j] = '0';
//         }
//         code[dim][j] = '0' + (code[dim][j] - '0' + 1) % 3;
//     }
// }


// template<std::size_t N>
// void NSequential::Fragment<N>::DimInvTransformToPointCode(std::size_t dim, bool leftBorder) {
//     if (
//         (leftBorder && isOddCentral[dim]) || (!leftBorder && !isOddCentral[dim])
//     ) {
//         std::size_t j = code[dim].size() - 1;
//         for (; j != 0 && code[dim][j] == '0'; --j) {
//             code[dim][j] = '2';
//         }
//         code[dim][j] = '0' + (code[dim][j] - '0' + 2) % 3;
//     }
// }


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
void NSequential::Fragment<N>::traversePointTree(char planeId) {
    // std::cout << "planeId = " << planeId << '\n';
    std::shared_ptr<NSequential::OneDimPointTree> from;
    std::size_t zero_count;
    if (isOddCentral[prevDivDim]) {
        // std::cout << "from right\n";
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
        // std::cout << "from left\n";
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
    // std::cout << "zero_cout = " << zero_count << '\n';
    // std::cout << "X_left = ";
    // for (auto el : leftPointTreeRefs) {
    //     std::cout << el->x << ' ';
    // }
    // std::cout << '\n';
    // std::cout << "X_right = ";
    // for (auto el : rightPointTreeRefs) {
    //     std::cout << el->x << ' ';
    // }
    // std::cout << '\n';
    // std::cout << "right = " << rightZeroDivsCount[prevDivDim] << ' ';
    // std::cout << "left = " << leftZeroDivsCount[prevDivDim] << '\n';

}


template<std::size_t N>
void NSequential::Fragment<N>::Divide(char planeId) {
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
    // std::cout << "New code = " << code[divDim] << '\n' << std::endl;

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


