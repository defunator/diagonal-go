#pragma once
#include "fragment.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include <chrono>


namespace NSequential {

template<std::size_t N>
class Plane {
private:
    // but this structure may be not memory efficient
    // because there could be 2 ^ N records for each point
    // May be better use prefix tree or int128_t (but this way has also limits with precision eps <= 1 / 128)
    // or K-d tree but this will have costly operations with double
    std::function<double(const std::array<double, N>&)> f;
    std::array<double, N> left;
    std::array<double, N> right;
    std::vector<std::array<double, N>> pointCoordinates;
    std::vector<double> pointValues;
    double lambdaMax;
    std::size_t bestPoint;
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
    std::map<std::string, std::multiset<std::size_t>> codeToPointId[N];
    std::vector<NSequential::Fragment<N>> searchFragments; 

    Plane(
        const std::array<double, N>& left,
        const std::array<double, N>& right,
        const std::function<double(const std::array<double, N>&)>& f,
        double r = 1.1,
        double C = 100
    );

    std::size_t GetPointIdByCode(const std::array<std::string, N> &pointCode);
    std::size_t GetPointIdByFragmentCode(std::size_t fragmentId, bool leftBorder = true);

    void DivideFragment(std::size_t fragmentId);

    double GetBestFragmentDiff();
    double GetBestPoint(std::array<double, N>& optimum);

    std::size_t GetBestFragmentId();

    std::size_t FCount();
    
    bool AssertFragment(std::size_t fragmentId) {
        bool fail = false;
        double eps = 1e-4;
        searchFragments[fragmentId].TransformToPointCode();
        std::array<double, N> x;
        for (std::size_t i = 0; i != N; ++i) {
            x[i] = left[i];
            double p = 1.;
            for (auto c : searchFragments[fragmentId].code[i]) {
                double w = (c - '0');
                x[i] += (right[i] - left[i]) * w / p;
                p *= 3.;
            }
        }
        searchFragments[fragmentId].InvTransformToPointCode();
        if (std::abs(f(x) - pointValues[searchFragments[fragmentId].leftPointId]) > eps) {
            std::cout << "++++++++++++++++++++\n";
            std::cout << "fragment code = ";
            for (auto code : searchFragments[fragmentId].code) {
                std::cout << "(" << code << ") ";
            }
            std::cout << '\n';
            searchFragments[fragmentId].TransformToPointCode();
            std::cout << "left point code = ";
            for (auto code : searchFragments[fragmentId].code) {
                std::cout << "(" << code << ") ";
            }
            std::cout << '\n';
            std::cout << "f( ";
            for (auto el : x) {
                std::cout << el << ' ';
            }
            std::cout << ") = " << f(x) << " != " << pointValues[searchFragments[fragmentId].leftPointId] << " = f( ";
            for (auto el : pointCoordinates[searchFragments[fragmentId].leftPointId]) {
                std::cout << el << ' ';
            }
            std::cout << ")\n++++++++++++++++++++\n";
            searchFragments[fragmentId].InvTransformToPointCode();
            fail = true;
        }

        searchFragments[fragmentId].TransformToPointCode(false);
        for (std::size_t i = 0; i != N; ++i) {
            x[i] = left[i];
            double p = 1.;
            for (auto c : searchFragments[fragmentId].code[i]) {
                double w = (c - '0');
                x[i] += (right[i] - left[i]) * w / p;
                p *= 3.;
            }
        }
        searchFragments[fragmentId].InvTransformToPointCode(false);
        if (std::abs(f(x) - pointValues[searchFragments[fragmentId].rightPointId]) > eps) {
            std::cout << "++++++++++++++++++++\n";
            std::cout << "fragment code = ";
            for (auto code : searchFragments[fragmentId].code) {
                std::cout << "(" << code << ") ";
            }
            std::cout << '\n';
            searchFragments[fragmentId].TransformToPointCode(false);
            std::cout << "right point code = ";
            for (auto code : searchFragments[fragmentId].code) {
                std::cout << "(" << code << ") ";
            }
            std::cout << '\n';
            std::cout << "f( ";
            for (auto el : x) {
                std::cout << el << ' ';
            }
            std::cout << ") = " << f(x) << " != " << pointValues[searchFragments[fragmentId].rightPointId] << " = f( ";
            for (auto el : pointCoordinates[searchFragments[fragmentId].rightPointId]) {
                std::cout << el << ' ';
            }
            std::cout << ")\n++++++++++++++++++++\n";
            searchFragments[fragmentId].InvTransformToPointCode(false);
            fail = true;
        }
        return fail;
    }
};

} // end namespace NSequential


template <std::size_t N>
NSequential::Plane<N>::Plane(
    const std::array<double, N> &left,
    const std::array<double, N> &right,
    const std::function<double(const std::array<double, N>&)> &f,
    double r,
    double C
) : f(f)
    , left(left)
    , right(right)
    , bestPoint(0)
    , r(r)
    , C(C)
{
    pointCoordinates.push_back(left);
    pointCoordinates.push_back(right);
    pointValues.push_back(f(left));
    pointValues.push_back(f(right));
    if (pointValues[1] < pointValues[0]) {
        ++bestPoint;
    }

    std::array<std::string, N> code;
    code.fill("1");
    addPoint(code, 1);
    code.fill("0"); 
    addPoint(code, 0);
    double diff = 0;
    std::priority_queue<std::pair<double, std::size_t>> diffDim;
    for (std::size_t i = 0; i != N; ++i) {
        diffDim.push({right[i] - left[i], i});
        assert(right[i] - left[i] > 0);
        diff += (right[i] - left[i]) * (right[i] - left[i]);
    }
    diff = std::sqrt(diff);
    searchFragments.emplace_back(std::move(code), std::move(diffDim), diff);
    searchFragments.back().leftPointId = 0, searchFragments.back().rightPointId = 1;
    lambdaMax = std::max(
        std::abs(pointValues[0] - pointValues[1]) / searchFragments[0].diff,
        std::numeric_limits<double>::epsilon()
    );
}


template <std::size_t N>
std::size_t NSequential::Plane<N>::getLastNonZero(const std::string code) {
    std::size_t res = code.size() - 1;
    for (; res != 0 && code[res] == '0'; --res);
    return res;
}


template <std::size_t N>
void NSequential::Plane<N>::addPoint(
    const std::array<std::string, N>& pointCode,
    std::size_t pointId
) {
    for (std::size_t i = 0; i != N; ++i) {
        std::size_t lastNonZero = getLastNonZero(pointCode[i]);
        std::string realCode = pointCode[i].substr(0, lastNonZero + 1);
        codeToPointId[i][realCode].insert(pointId);
    }
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
    deletePoint(searchFragments[fragmentId].code);
    searchFragments[fragmentId].InvTransformToPointCode();
    searchFragments[fragmentId].TransformToPointCode(false);
    deletePoint(searchFragments[fragmentId].code);
    searchFragments[fragmentId].InvTransformToPointCode(false);
}


template<std::size_t N>
void NSequential::Plane<N>::addFragmentPoints(
    std::size_t fragmentId,
    std::size_t left,
    std::size_t right
) {
    if (left != 0) {
        searchFragments[fragmentId].TransformToPointCode();
        addPoint(searchFragments[fragmentId].code, left-1);
        searchFragments[fragmentId].InvTransformToPointCode();
    }
    if (right != 0) {
        searchFragments[fragmentId].TransformToPointCode(false);
        addPoint(searchFragments[fragmentId].code, right-1);
        searchFragments[fragmentId].InvTransformToPointCode(false);
    }
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
    std::size_t pointId = GetPointIdByCode(searchFragments[fragmentId].code);
    searchFragments[fragmentId].InvTransformToPointCode(leftBorder);
    return pointId;
}


template<std::size_t N>
void NSequential::Plane<N>::DivideFragment(std::size_t fragmentId) {
    std::size_t leftBase = searchFragments[fragmentId].leftPointId;
    std::size_t rightBase = searchFragments[fragmentId].rightPointId;

    searchFragments.push_back(searchFragments[fragmentId]);
    searchFragments.back().Divide('0');
    searchFragments.push_back(searchFragments[fragmentId]);
    searchFragments.back().Divide('2');
    searchFragments[fragmentId].Divide('1');

    std::size_t left = GetPointIdByFragmentCode(fragmentId);
    std::size_t right = GetPointIdByFragmentCode(fragmentId, false);

    if (left == 0) {
        searchFragments[fragmentId].TransformToPointCode();
        std::array<double, N> leftCoord;
        for (std::size_t i = 0; i != N; ++i) {
            leftCoord[i] = pointCoordinates[0][i];
            double p = 1.;
            for (auto c : searchFragments[fragmentId].code[i]) {
                double w = (c - '0');
                leftCoord[i] += (pointCoordinates[1][i] - pointCoordinates[0][i]) * w / p;
                p *= 3.;
            }
        }
        pointValues.push_back(f(leftCoord));
        if (pointValues.back() < pointValues[bestPoint]) {
            bestPoint = pointValues.size() - 1;
        }
        searchFragments[fragmentId].InvTransformToPointCode();

        pointCoordinates.push_back(std::move(leftCoord));
        left = pointCoordinates.size();
    }

    if (right == 0) {
        searchFragments[fragmentId].TransformToPointCode(false);
        std::array<double, N> rightCoord;
        for (std::size_t i = 0; i != N; ++i) {
            rightCoord[i] = pointCoordinates[0][i];
            double p = 1.;
            for (auto c : searchFragments[fragmentId].code[i]) {
                double w = (c - '0');
                rightCoord[i] += (pointCoordinates[1][i] - pointCoordinates[0][i]) * w / p;
                p *= 3.;
            }
        }
        searchFragments[fragmentId].InvTransformToPointCode(false);

        pointValues.push_back(f(rightCoord));
        if (pointValues.back() < pointValues[bestPoint]) {
            bestPoint = pointValues.size() - 1;
        }
        pointCoordinates.push_back(std::move(rightCoord));
        right = pointCoordinates.size();
    }

    addFragmentPoints(fragmentId, left, right);
    addFragmentPoints(searchFragments.size() - 1, left, 0);
    addFragmentPoints(searchFragments.size() - 2, 0, right);

    --left, --right;
    std::vector<std::size_t> pointsToTest({left, right, rightBase, leftBase});
    for (auto x1 : pointsToTest) {
        for (auto x2 : pointsToTest) {
            lambdaMax = std::max(
                lambdaMax,
                std::abs(pointValues[x1] - pointValues[x2]) / searchFragments.back().diff
            );
        }
    }
    searchFragments[fragmentId].leftPointId = left, searchFragments[fragmentId].rightPointId = right;
    searchFragments[searchFragments.size() - 1].leftPointId = left, searchFragments[searchFragments.size() - 1].rightPointId = rightBase;
    searchFragments[searchFragments.size() - 2].leftPointId = leftBase, searchFragments[searchFragments.size() - 2].rightPointId = right;
}


template<std::size_t N>
double NSequential::Plane<N>::GetBestFragmentDiff() {
    std::size_t bestI = 0;
    for (std::size_t i = 0; i != searchFragments.size(); ++i) {
        if (searchFragments[bestI].R < searchFragments[i].R) {
            bestI = i;
        }
    }
    return searchFragments[bestI].diff;
}


template<std::size_t N>
double NSequential::Plane<N>::GetBestPoint(std::array<double, N>& optimum) {
    optimum = pointCoordinates[bestPoint];
    return pointValues[bestPoint];
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::GetBestFragmentId() {
    std::size_t k = std::max(searchFragments.size() / 2, 1ul);
    double mu = (r + C / k) * lambdaMax;
    std::size_t bestFragment = 0;
    for (std::size_t i = 0; i != searchFragments.size(); ++i) {
        // if (AssertFragment(i)) {
        //     std::cout << searchFragments[i].leftPointId << ' ' << searchFragments[i].rightPointId << std::endl;
        //     exit(0);
        // }
        double fLeft = pointValues[searchFragments[i].leftPointId];
        double fRight = pointValues[searchFragments[i].rightPointId];
        searchFragments[i].UpdR(mu, fLeft, fRight);
        if (searchFragments[bestFragment].R < searchFragments[i].R) {
            bestFragment = i;
        }
    }
    return bestFragment;
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::FCount() {
    return pointValues.size();
}

