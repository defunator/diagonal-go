// tree plane
#pragma once
#include "fragment.hpp"
// #include "point_tree.hpp"

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
    PointTree<N> codeToPointId;
    double r;
    double C;
    double lambdaMax;
    std::size_t bestPoint;
    std::size_t bestFragmentId;

    // void addFragmentPointCodes(std::size_t fragmentId, std::size_t left, std::size_t right);
    void tryAddFragmentPoints(std::size_t fragmentId);

public:
    std::vector<NSequential::Fragment<N>> searchFragments; 

    Plane(
        const std::array<double, N>& left,
        const std::array<double, N>& right,
        const std::function<double(const std::array<double, N>&)>& f,
        double r = 1.1,
        double C = 100
    );

    void DivideFragment(std::size_t fragmentId);

    double GetBestFragmentDiff();
    double GetBestPoint(std::array<double, N>& optimum);

    std::size_t GetBestFragmentId();

    std::size_t FCount();

    // bool AssertFragment(std::size_t fragmentId) {
    //     bool fail = false;
    //     double eps = 1e-4;
    //     searchFragments[fragmentId].TransformToPointCode();
    //     std::array<double, N> x;
    //     for (std::size_t i = 0; i != N; ++i) {
    //         x[i] = codeToPointId.left[i];
    //         double p = 1.;
    //         for (auto c : searchFragments[fragmentId].code[i]) {
    //             double w = (c - '0');
    //             x[i] += (codeToPointId.right[i] - codeToPointId.left[i]) * w / p;
    //             p *= 3.;
    //         }
    //     }
    //     searchFragments[fragmentId].InvTransformToPointCode();
    //     double diff = 0;
    //     for (std::size_t i = 0; i != N; ++i) {
    //         diff += (x[i] - codeToPointId.points[searchFragments[fragmentId].leftPointId-1][i]) * (x[i] - codeToPointId.points[searchFragments[fragmentId].leftPointId-1][i]);
    //     }
    //     if (std::abs(codeToPointId.f(x) - searchFragments[fragmentId].fLeft) > eps || diff > eps) {
    //         std::cout << "++++++++++++++++++++\n";
    //         std::cout << "fragment code = ";
    //         for (auto code : searchFragments[fragmentId].code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         searchFragments[fragmentId].TransformToPointCode();
    //         std::cout << "left point code = ";
    //         for (auto code : searchFragments[fragmentId].code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         std::cout << "f( ";
    //         for (auto el : x) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ") = " << codeToPointId.f(x) << " != " << searchFragments[fragmentId].fLeft;
    //         std::cout << " = f( ";
    //         for (auto el : codeToPointId.points[searchFragments[fragmentId].leftPointId-1]) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ")\n";
    //         std::cout << "++++++++++++++++++++\n";
    //         searchFragments[fragmentId].InvTransformToPointCode();
    //         fail = true;
    //     }

    //     searchFragments[fragmentId].TransformToPointCode(false);
    //     for (std::size_t i = 0; i != N; ++i) {
    //         x[i] = codeToPointId.left[i];
    //         double p = 1.;
    //         for (auto c : searchFragments[fragmentId].code[i]) {
    //             double w = (c - '0');
    //             x[i] += (codeToPointId.right[i] - codeToPointId.left[i]) * w / p;
    //             p *= 3.;
    //         }
    //     }
    //     searchFragments[fragmentId].InvTransformToPointCode(false);
    //     diff = 0;
    //     for (std::size_t i = 0; i != N; ++i) {
    //         diff += (x[i] - codeToPointId.points[searchFragments[fragmentId].rightPointId-1][i]) * (x[i] - codeToPointId.points[searchFragments[fragmentId].rightPointId-1][i]);
    //     }
    //     if (std::abs(codeToPointId.f(x) - searchFragments[fragmentId].fRight) > eps || diff > eps) {
    //         std::cout << "++++++++++++++++++++\n";
    //         std::cout << "fragment code = ";
    //         for (auto code : searchFragments[fragmentId].code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         searchFragments[fragmentId].TransformToPointCode(false);
    //         std::cout << "right point code = ";
    //         for (auto code : searchFragments[fragmentId].code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         std::cout << "f( ";
    //         for (auto el : x) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ") = " << codeToPointId.f(x) << " != " << searchFragments[fragmentId].fRight;
    //         std::cout << " = f( ";
    //         for (auto el : codeToPointId.points[searchFragments[fragmentId].rightPointId-1]) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ")\n";
    //         std::cout << "++++++++++++++++++++\n";
    //         searchFragments[fragmentId].InvTransformToPointCode(false);
    //         fail = true;
    //     }
    //     return fail;
    // }
};

} // end namespace NSequential


template <std::size_t N>
NSequential::Plane<N>::Plane(
    const std::array<double, N>& left,
    const std::array<double, N>& right,
    const std::function<double(const std::array<double, N>&)>& f,
    double r,
    double C
)
    : codeToPointId(left, right, f)
    , r(r)
    , C(C)
{
    std::array<std::string, N> code;
    code.fill("0");
    std::array<std::shared_ptr<OneDimPointTree>, N> leftPointTreeRefs;
    codeToPointId.GetTreeNodesByCode(code, leftPointTreeRefs);
    code.fill("1");
    std::array<std::shared_ptr<OneDimPointTree>, N> rightPointTreeRefs;
    codeToPointId.GetTreeNodesByCode(code, rightPointTreeRefs);

    searchFragments.reserve(50000);
    searchFragments.emplace_back(
        left, right,
        std::move(leftPointTreeRefs), std::move(rightPointTreeRefs),
        codeToPointId.GetValue(1), codeToPointId.GetValue(2)
    );
    lambdaMax = std::max(
        std::abs(codeToPointId.GetValue(1) - codeToPointId.GetValue(2)) / searchFragments[0].diff,
        std::numeric_limits<double>::epsilon()
    );
}


// template<std::size_t N>
// void NSequential::Plane<N>::addFragmentPointCodes(
//     std::size_t fragmentId,
//     std::size_t left,
//     std::size_t right
// ) {
//     searchFragments[fragmentId].TransformToPointCode();
//     codeToPointId.AddPointCode(searchFragments[fragmentId].code, left);
//     searchFragments[fragmentId].InvTransformToPointCode();

//     searchFragments[fragmentId].TransformToPointCode(false);
//     codeToPointId.AddPointCode(searchFragments[fragmentId].code, right);
//     searchFragments[fragmentId].InvTransformToPointCode(false);
// }


template <std::size_t N>
void NSequential::Plane<N>::tryAddFragmentPoints(std::size_t fragmentId) {
    std::size_t divDim = searchFragments[fragmentId].prevDivDim;
    // searchFragments[fragmentId].DimTransformToPointCode(divDim);
    // codeToPointId.GetTreeNodeByCode(
    //     searchFragments[fragmentId].code[divDim],
    //     searchFragments[fragmentId].leftPointTreeRefs[divDim],
    //     divDim
    // );
    std::size_t left = codeToPointId.GetPointIdFromArray(searchFragments[fragmentId].leftPointTreeRefs);
    if (left == 0) {
        left = codeToPointId.AddPoint(
            searchFragments[fragmentId].leftPointTreeRefs
        );
    }
    searchFragments[fragmentId].fLeft = codeToPointId.GetValue(left);
    // searchFragments[fragmentId].DimInvTransformToPointCode(divDim);

    // searchFragments[fragmentId].DimTransformToPointCode(divDim, false);
    // codeToPointId.GetTreeNodeByCode(
    //     searchFragments[fragmentId].code[divDim],
    //     searchFragments[fragmentId].rightPointTreeRefs[divDim],
    //     divDim
    // );
    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::size_t right = codeToPointId.GetPointIdFromArray(searchFragments[fragmentId].rightPointTreeRefs);
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // std::cout << "TIME = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;
    if (right == 0) {
        right = codeToPointId.AddPoint(
            searchFragments[fragmentId].rightPointTreeRefs
        );
    }
    searchFragments[fragmentId].fRight = codeToPointId.GetValue(right);
    // searchFragments[fragmentId].DimInvTransformToPointCode(divDim, false);

    searchFragments[fragmentId].leftPointId = left;
    searchFragments[fragmentId].rightPointId = right;
    searchFragments[searchFragments.size() - 1].leftPointId = left;
    searchFragments[searchFragments.size() - 2].rightPointId = right;

    searchFragments[searchFragments.size() - 1].fLeft = searchFragments[fragmentId].fLeft;
    searchFragments[searchFragments.size() - 2].fRight = searchFragments[fragmentId].fRight;
    searchFragments[searchFragments.size() - 1].leftPointTreeRefs[divDim] = searchFragments[fragmentId].leftPointTreeRefs[divDim];
    searchFragments[searchFragments.size() - 2].rightPointTreeRefs[divDim] = searchFragments[fragmentId].rightPointTreeRefs[divDim];

    codeToPointId.AddPointCode(searchFragments[fragmentId].leftPointTreeRefs, left);
    codeToPointId.AddPointCode(searchFragments[fragmentId].rightPointTreeRefs, right);
    codeToPointId.AddPointCode(searchFragments[searchFragments.size() - 1].leftPointTreeRefs, left);
    codeToPointId.AddPointCode(searchFragments[searchFragments.size() - 2].rightPointTreeRefs, right);
}


template<std::size_t N>
void NSequential::Plane<N>::DivideFragment(std::size_t fragmentId) {
    searchFragments.push_back(searchFragments[fragmentId]);
    searchFragments.back().Divide('0');
    searchFragments.push_back(searchFragments[fragmentId]);
    searchFragments.back().Divide('2');
    searchFragments[fragmentId].Divide('1');
    
    tryAddFragmentPoints(fragmentId);
    
    double f1 = searchFragments[searchFragments.size() - 2].fRight;
    double f2 = searchFragments[searchFragments.size() - 2].fLeft;
    lambdaMax = std::max(lambdaMax, std::abs(f1 - f2) / searchFragments.back().diff);

    f1 = searchFragments[searchFragments.size() - 1].fRight;
    f2 = searchFragments[searchFragments.size() - 1].fLeft;
    lambdaMax = std::max(lambdaMax, std::abs(f1 - f2) / searchFragments.back().diff);

    f1 = searchFragments[fragmentId].fRight;
    f2 = searchFragments[fragmentId].fLeft;
    lambdaMax = std::max(lambdaMax, std::abs(f1 - f2) / searchFragments.back().diff);

    double diff = 2. * searchFragments.back().prevDiff / 3.;
    f1 = searchFragments[searchFragments.size() - 2].fLeft;
    f2 = searchFragments[fragmentId].fLeft;
    lambdaMax = std::max(lambdaMax, std::abs(f1 - f2) / diff);

    f1 = searchFragments[searchFragments.size() - 1].fRight;
    f2 = searchFragments[fragmentId].fRight;
    lambdaMax = std::max(lambdaMax, std::abs(f1 - f2) / diff);


}


template<std::size_t N>
double NSequential::Plane<N>::GetBestFragmentDiff() {
    return searchFragments[bestFragmentId].diff;
}


template<std::size_t N>
double NSequential::Plane<N>::GetBestPoint(std::array<double, N>& optimum) {
    return codeToPointId.GetBestPoint(optimum);
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::GetBestFragmentId() {
    std::size_t k = std::max(searchFragments.size() / 2, 1ul);
    double mu = (r + C / k) * lambdaMax;
    std::size_t bestFragment = 0;
    for (std::size_t i = 0; i != searchFragments.size(); ++i) {
        // if (AssertFragment(i)) {
        //     exit(0);
        // }
        double fLeft = searchFragments[i].fLeft;
        double fRight = searchFragments[i].fRight;
        if (fLeft > fRight) {
            std::swap(fLeft, fRight);
        }
        searchFragments[i].UpdR(mu, fLeft, fRight);
        if (searchFragments[bestFragment].R < searchFragments[i].R) {
            bestFragment = i;
        }
    }
    bestFragmentId = bestFragment;
    return bestFragment;
}


template<std::size_t N>
std::size_t NSequential::Plane<N>::FCount() {
    return codeToPointId.FCount();
}

