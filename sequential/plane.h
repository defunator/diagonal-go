#pragma once
#include "fragment.h"

#include <array>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
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

}


#include "plane.tpp"
