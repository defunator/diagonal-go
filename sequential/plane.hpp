#pragma once
#include "fragment.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>  // std::numeric_limits<double>::epsilon()
#include <functional>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


namespace NSequential {

template<std::size_t N>
class Plane {
private:
    // but this structure may be not memory efficient
    // because there could be 2 ^ N records for each point
    // May be better use prefix tree or int128_t (but this way has also limits with precision eps <= 1 / 128)
    // or K-d tree but this will have costly operations with double
    std::unordered_map<std::string, std::multiset<std::size_t>> codeToPointId[N];
    std::vector<NSequential::Fragment<N>> searchFragments; 
    std::vector<std::array<double, N>> pointCoordinates;
    std::vector<double> pointValues;
    std::function<double(std::array<double, N>)> f;
    double lambdaMax;
    std::size_t bestPoint;
    double eps;
    double r;
    double C;

    std::size_t getLastNonZero(const std::string code);

    void addPoint(const std::array<std::string, N>& pointCode, std::size_t pointId);

    void deletePoint(const std::array<std::string, N>& pointCode);

    void updPointIds(std::multiset<std::size_t> &pointIds,
                    const std::multiset<std::size_t> &newPointIds);

    void deleteFragmentPoints(std::size_t fragmentId);

    void addFragmentPoints(std::size_t fragmentId, std::size_t left, std::size_t right);

public:
    Plane(
        const std::array<double, N>& left,
        const std::array<double, N>& right,
        const std::function<double(std::array<double, N>)>& f,
        double eps = 0.01,
        double r = 1.1,
        double C = 100
    );

    std::size_t GetPointIdByCode(const std::array<std::string, N> &pointCode);
    std::size_t GetPointIdByFragmentCode(std::size_t fragmentId, bool leftBorder = true);

    void DivideFragment(std::size_t fragmentId);

    double GetBestFragmentDiff2();
    double GetBestPoint(std::array<double, N>& optimum);

    std::size_t GetBestFragmentId();
};

} // end namespace NSequential


template<std::size_t N>
std::size_t NSequential::Plane<N>::getLastNonZero(const std::string code) {
    std::size_t res = 0;
    for (std::size_t i = 0; i != code.size(); ++i) {
        if (code[i] != '0') {
            res = i;
        }
    }
    return res;
}


template<std::size_t N>
void NSequential::Plane<N>::addPoint(
    const std::array<std::string, N>& pointCode,
    std::size_t pointId
) {
    for (std::size_t i = 0; i != N; ++i) {
        std::size_t lastNonZero = getLastNonZero(pointCode[i]);
        std::string realCode = pointCode[i].substr(0, lastNonZero + 1);
        codeToPointId[i][realCode].insert(pointId);
    }
    assert(GetPointIdByCode(pointCode) != 0);
}


template<std::size_t N>
void NSequential::Plane<N>::deletePoint(
    const std::array<std::string, N>& pointCode
) {
    std::size_t pointId = GetPointIdByCode(pointCode) - 1;
    for (std::size_t i = 0; i != N; ++i) {
        const std::string &code = pointCode[i];
        std::size_t lastNonZero = getLastNonZero(code);
        std::string realCode = code.substr(0, lastNonZero + 1);
        auto idToRemove = codeToPointId[i][realCode].find(pointId);
        codeToPointId[i][realCode].erase(idToRemove);
        if (codeToPointId[i][realCode].empty()) {
            codeToPointId[i].erase(realCode);
        }
    }
}


template<std::size_t N>
void NSequential::Plane<N>::updPointIds(
    std::multiset<std::size_t> &pointIds,
    const std::multiset<std::size_t> &newPointIds
) {
    if (pointIds.empty()) {
        pointIds = newPointIds;
    } else {
        std::vector<std::size_t> tmp;
        std::set_intersection(
            pointIds.begin(), pointIds.end(),
            newPointIds.begin(), newPointIds.end(),
            std::back_inserter(tmp));
        pointIds.clear();
        for (auto& el : tmp) {
            pointIds.insert(el);
        }
    }
}


template<std::size_t N>
void NSequential::Plane<N>::deleteFragmentPoints(std::size_t fragmentId) {
    searchFragments[fragmentId].TransformToPointCode();
    // std::cout << "Try delete" << std::endl;
    deletePoint(searchFragments[fragmentId].getCode());
    searchFragments[fragmentId].InvTransformToPointCode();
    searchFragments[fragmentId].TransformToPointCode(false);
    deletePoint(searchFragments[fragmentId].getCode());
    searchFragments[fragmentId].InvTransformToPointCode(false);
}


template<std::size_t N>
void NSequential::Plane<N>::addFragmentPoints(
    std::size_t fragmentId,
    std::size_t left,
    std::size_t right
) {
    searchFragments[fragmentId].TransformToPointCode();
    addPoint(searchFragments[fragmentId].getCode(), left);
    searchFragments[fragmentId].InvTransformToPointCode();
    searchFragments[fragmentId].TransformToPointCode(false);
    addPoint(searchFragments[fragmentId].getCode(), right);
    searchFragments[fragmentId].InvTransformToPointCode(false);
}


template<std::size_t N>
NSequential::Plane<N>::Plane(
    const std::array<double, N> &left,
    const std::array<double, N> &right,
    const std::function<double(std::array<double, N>)> &f,
    double eps,
    double r,
    double C
) : f(f)
    , bestPoint(0)
    , eps(eps)
    , r(r)
    , C(C)
{
    pointCoordinates.push_back(left);
    pointValues.push_back(f(left)); // async?
    pointCoordinates.push_back(right);
    pointValues.push_back(f(right));
    if (pointValues[1] < pointValues[0]) {
        ++bestPoint;
    }
    std::array<std::string, N> code;
    code.fill("1");
    addPoint(code, 0);
    code.fill("0"); 
    addPoint(code, 1);
    std::set<std::pair<double, std::size_t>> diff;
    for (std::size_t i = 0; i != N; ++i) {
        diff.insert({right[i] - left[i], i});
        assert(right[i] - left[i] > 0);
    }
    searchFragments.push_back(std::move(
        NSequential::Fragment<N>(std::move(code), std::move(diff))
    ));
    lambdaMax = std::max(
        std::abs(pointValues[0] - pointValues[1]) / std::sqrt(searchFragments[0].diff2),
        std::numeric_limits<double>::epsilon()
    );
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::GetPointIdByCode(
    const std::array<std::string, N> &pointCode
) {
    std::multiset<std::size_t> pointIds;
    for (std::size_t i = 0; i != N; ++i) {
        std::size_t lastNonZero = getLastNonZero(pointCode[i]);
        std::string realCode = pointCode[i].substr(0, lastNonZero + 1);
        updPointIds(pointIds, codeToPointId[i][realCode]);
        if (pointIds.empty()) { // Point not found
            return 0;
        }
    }
    return *pointIds.begin() + 1;
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::GetPointIdByFragmentCode(
    std::size_t fragmentId,
    bool leftBorder
) {
    searchFragments[fragmentId].TransformToPointCode(leftBorder);
    std::size_t pointId = GetPointIdByCode(searchFragments[fragmentId].getCode());
    searchFragments[fragmentId].InvTransformToPointCode(leftBorder);
    return pointId;
}


template<std::size_t N>
void NSequential::Plane<N>::DivideFragment(std::size_t fragmentId) {
    std::size_t divideDim = searchFragments[fragmentId].getDivideDim();
    std::size_t leftBase = GetPointIdByFragmentCode(fragmentId) - 1;
    std::size_t rightBase = GetPointIdByFragmentCode(fragmentId, false) - 1;

    deleteFragmentPoints(fragmentId);
    NSequential::Fragment<N> zeroFrag(searchFragments[fragmentId]);
    searchFragments[fragmentId].Divide('0', zeroFrag);
    searchFragments.push_back(std::move(zeroFrag));
    NSequential::Fragment<N> twoFrag(searchFragments[fragmentId]);
    searchFragments[fragmentId].Divide('2', twoFrag);
    searchFragments.push_back(std::move(twoFrag));
    searchFragments[fragmentId].InplaceDivide('1');

    std::size_t left = GetPointIdByFragmentCode(fragmentId);
    std::size_t right = GetPointIdByFragmentCode(fragmentId, false);
    if (left == 0) {
        std::array<double, N> leftCoord = pointCoordinates[leftBase];
        if (searchFragments[fragmentId].isOddDim(divideDim)) {
            leftCoord[divideDim] = (pointCoordinates[leftBase][divideDim] +
                                    2 * pointCoordinates[rightBase][divideDim]) / 3;
        } else {
            leftCoord[divideDim] = (2 * pointCoordinates[leftBase][divideDim] +
                                    pointCoordinates[rightBase][divideDim]) / 3;
        }
        pointValues.push_back(f(leftCoord));
        if (pointValues.back() < pointValues[bestPoint]) {
            bestPoint = pointValues.size() - 1;
        }
        pointCoordinates.push_back(std::move(leftCoord));
        left = pointCoordinates.size();
    }
    if (right == 0) {
        std::array<double, N> rightCoord = pointCoordinates[rightBase];
        if (searchFragments[fragmentId].isOddDim(divideDim)) {
            rightCoord[divideDim] = (2 * pointCoordinates[leftBase][divideDim] +
                                    pointCoordinates[rightBase][divideDim]) / 3;
        } else {
            rightCoord[divideDim] = (pointCoordinates[leftBase][divideDim] +
                                    2 * pointCoordinates[rightBase][divideDim]) / 3;
        }
        pointValues.push_back(f(rightCoord));
        if (pointValues.back() < pointValues[bestPoint]) {
            bestPoint = pointValues.size() - 1;
        }
        pointCoordinates.push_back(std::move(rightCoord));
        right = pointCoordinates.size();
    }
    --left, --right;
    addFragmentPoints(fragmentId, left, right);
    addFragmentPoints(searchFragments.size() - 1, left, rightBase);
    addFragmentPoints(searchFragments.size() - 2, leftBase, right);

    std::vector<std::size_t> pointsToTest({left, right, rightBase, leftBase});
    for (auto x1 : pointsToTest) {
        for (auto x2 : pointsToTest) {
            if (x1 == x2) {
                continue;
            }
            lambdaMax = std::max(
                lambdaMax,
                std::abs(pointValues[x1] - pointValues[x2]) / std::sqrt(searchFragments.back().diff2)
            );
        }
    }
}


template<std::size_t N>
double NSequential::Plane<N>::GetBestFragmentDiff2() {
    std::size_t bestI = 0;
    for (std::size_t i = 0; i != searchFragments.size(); ++i) {
        if (searchFragments[bestI].R < searchFragments[i].R) {
            bestI = i;
        }
    }
    return searchFragments[bestI].diff2;
}


template<std::size_t N>
double NSequential::Plane<N>::GetBestPoint(std::array<double, N>& optimum) {
    optimum = pointCoordinates[bestPoint];
    return pointValues[bestPoint];
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::GetBestFragmentId() {
    // double dMax = 0;
    // double lambdaMax = 0;
    // for (std::size_t i = 0; i != searchFragments.size(); ++i) {
    //     dMax = std::max(dMax, searchFragments[i].diff2);
    //     double fLeft = pointValues[GetPointIdByFragmentCode(i)];
    //     double fRight = pointValues[GetPointIdByFragmentCode(i, false)];
    //     lambdaMax = std::max(
    //         lambdaMax,
    //         std::abs(fRight - fLeft) / std::sqrt(searchFragments[i].diff2)
    //     );
    // }
    std::cout << "lambdaMax = " << lambdaMax << std::endl;
    double dMax = std::sqrt(dMax);
    std::size_t bestFragment = 0;
    std::size_t k = searchFragments.size() / 2;
    for (std::size_t i = 0; i != searchFragments.size(); ++i) {
        double fLeft = pointValues[GetPointIdByFragmentCode(i)];
        double fRight = pointValues[GetPointIdByFragmentCode(i, false)];
        searchFragments[i].updR(C, r, k, dMax, lambdaMax, fLeft, fRight);
        if (searchFragments[bestFragment].R < searchFragments[i].R) {
            bestFragment = i;
        }
    }
    std::cout << "bestFragment = " << bestFragment << std::endl;
    std::cout << "diff = " << searchFragments[bestFragment].diff2 << std::endl;
    std::cout << "sz = " << pointValues.size() << std::endl;
    std::cout << "fragSZ = " << searchFragments.size() << std::endl;
    return bestFragment;
}
