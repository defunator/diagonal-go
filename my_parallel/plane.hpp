// tree plane
#pragma once
#include "atomic_list.hpp"
#include "fragment.hpp"
// #include "point_tree.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <tbb/concurrent_priority_queue.h>
#include <unordered_map>
#include <utility>
#include <vector>

#include <chrono>


namespace NParallel {

template <std::size_t N>
class Plane {
public:
    typedef std::shared_ptr<Node<Fragment<N>>> NodePtr;

private:
    PointTree<N> codeToPointId;
    double r;
    double C;
    std::atomic<double> lambdaMax;
    std::atomic<double> prevMu;
    // NodePtr bestFragment;
    tbb::concurrent_priority_queue<std::pair<double, NodePtr>> fragmentQueue;
    AtomicList<Fragment<N>> searchFragments;

    /**
     * Processes new three fragments (if needed eval new points etc)
     * You should ensure on your own that after `fragment` Node are another two fragments
     * from this split
     * 
     * @param fragment pointer in `searchFragments` linked list to first fragment after triple division
     */
    void tryAddFragmentPoints(NodePtr& fragment);

    void updateLambdaMax(double f1, double f2, double diff);

public:
    Plane(
        const std::array<double, N>& left,
        const std::array<double, N>& right,
        const std::function<double(const std::array<double, N>&)>& f,
        double r = 1.1,
        double C = 100
    );

    void DivideFragment(NodePtr& fragment);

    double GetBestPoint(std::array<double, N>& optimum);

    bool GetBestFragment(NodePtr& to);
    void PutFragmentBack(NodePtr& fragment);
    void RecalcAfterDivision(NodePtr& fragment);
    // Warning: swaps `fragmentQueue` and iterates `searchFragments`, so it is not thread-safe
    void RecalcAllFragments();

    std::size_t FCount() const { return codeToPointId.FCount(); }
    std::size_t NQueued() const { return fragmentQueue.size(); }
    std::size_t NFragments() const { return searchFragments.Size(); }

    double PrevUpdMuDiff() const {
        std::size_t k = std::max(searchFragments.Size() / 2, 1ul);
        double mu = (r + C / k) * lambdaMax;
        return std::abs(mu - prevMu);
    }

    // bool AssertFragment(Fragment<N>& fragment) {
    //     bool fail = false;
    //     double eps = 1e-4;
    //     fragment.TransformToPointCode();
    //     std::array<double, N> x;
    //     for (std::size_t i = 0; i != N; ++i) {
    //         x[i] = codeToPointId.left[i];
    //         double p = 1.;
    //         for (auto c : fragment.code[i]) {
    //             double w = (c - '0');
    //             x[i] += (codeToPointId.right[i] - codeToPointId.left[i]) * w / p;
    //             p *= 3.;
    //         }
    //     }
    //     fragment.InvTransformToPointCode();
    //     double diff = 0;
    //     // for (std::size_t i = 0; i != N; ++i) {
    //     //     diff += (x[i] - codeToPointId.points[fragment.leftPointId-1][i]) * (x[i] - codeToPointId.points[searchFragments[fragmentId].leftPointId-1][i]);
    //     // }
    //     if (std::abs(codeToPointId.f(x) - fragment.fLeft) > eps || diff > eps) {
    //         std::cout << "++++++++++++++++++++\n";
    //         std::cout << "fragment code = ";
    //         for (auto code : fragment.code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         fragment.TransformToPointCode();
    //         std::cout << "left point code = ";
    //         for (auto code : fragment.code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         std::cout << "f( ";
    //         for (auto el : x) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ") = " << codeToPointId.f(x) << " != " << fragment.fLeft;
    //         std::cout << " = f( ";
    //         // for (auto el : codeToPointId.points[searchFragments[fragmentId].leftPointId-1]) {
    //         //     std::cout << el << ' ';
    //         // }
    //         std::cout << ")\n";
    //         std::cout << "++++++++++++++++++++\n";
    //         fragment.InvTransformToPointCode();
    //         fail = true;
    //     }

    //     fragment.TransformToPointCode(false);
    //     for (std::size_t i = 0; i != N; ++i) {
    //         x[i] = codeToPointId.left[i];
    //         double p = 1.;
    //         for (auto c : fragment.code[i]) {
    //             double w = (c - '0');
    //             x[i] += (codeToPointId.right[i] - codeToPointId.left[i]) * w / p;
    //             p *= 3.;
    //         }
    //     }
    //     fragment.InvTransformToPointCode(false);
    //     diff = 0;
    //     // for (std::size_t i = 0; i != N; ++i) {
    //     //     diff += (x[i] - codeToPointId.points[searchFragments[fragmentId].rightPointId-1][i]) * (x[i] - codeToPointId.points[searchFragments[fragmentId].rightPointId-1][i]);
    //     // }
    //     if (std::abs(codeToPointId.f(x) - fragment.fRight) > eps || diff > eps) {
    //         std::cout << "++++++++++++++++++++\n";
    //         std::cout << "fragment code = ";
    //         for (auto code : fragment.code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         fragment.TransformToPointCode(false);
    //         std::cout << "right point code = ";
    //         for (auto code : fragment.code) {
    //             std::cout << "(" << code << ") ";
    //         }
    //         std::cout << '\n';
    //         std::cout << "f( ";
    //         for (auto el : x) {
    //             std::cout << el << ' ';
    //         }
    //         std::cout << ") = " << codeToPointId.f(x) << " != " << fragment.fRight;
    //         std::cout << " = f( ";
    //         // for (auto el : codeToPointId.points[fragment.rightPointId-1]) {
    //         //     std::cout << el << ' ';
    //         // }
    //         std::cout << ")\n";
    //         std::cout << "++++++++++++++++++++\n";
    //         fragment.InvTransformToPointCode(false);
    //         fail = true;
    //     }
    //     return fail;
    // }
};

} // end namespace NParallel


template <std::size_t N>
NParallel::Plane<N>::Plane(
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

    double fLeft = codeToPointId.GetPointIdFromArray(leftPointTreeRefs).second;
    double fRight = codeToPointId.GetPointIdFromArray(rightPointTreeRefs).second;
    searchFragments.InsertBack(std::move(Fragment<N>(
        left, right,
        std::move(leftPointTreeRefs), std::move(rightPointTreeRefs),
        fLeft, fRight
    )));
    NParallel::Plane<N>::NodePtr begin;
    searchFragments.GetBegin(begin);
    lambdaMax = std::max(
        std::abs(begin->value->fLeft - begin->value->fRight) / begin->value->diff,
        std::numeric_limits<double>::epsilon()
    );
    prevMu = lambdaMax * (r + C / 1.);
    RecalcAllFragments();
}


template <std::size_t N>
void NParallel::Plane<N>::tryAddFragmentPoints(
    NParallel::Plane<N>::NodePtr& fragment
) {
    std::shared_ptr<NParallel::Fragment<N>> oneFragment{fragment->value};
    std::shared_ptr<NParallel::Fragment<N>> zeroFragment{fragment->next->value};
    std::shared_ptr<NParallel::Fragment<N>> twoFragment{fragment->next->next->value};
    std::size_t divDim = oneFragment->prevDivDim;

    auto[leftId, leftValue] = codeToPointId.GetPointIdFromArray(oneFragment->leftPointTreeRefs);
    if (leftId == 0) {
        auto[newLeftId, newLeftValue] = codeToPointId.AddPoint(
            oneFragment->leftPointTreeRefs
        );
        std::swap(leftValue, newLeftValue);
    }
    oneFragment->fLeft = leftValue;

    auto[rightId, rightValue] = codeToPointId.GetPointIdFromArray(oneFragment->rightPointTreeRefs);
    if (rightId == 0) {
        auto[newRightId, newRightValue] = codeToPointId.AddPoint(
            oneFragment->rightPointTreeRefs
        );
        std::swap(rightValue, newRightValue);
    }
    oneFragment->fRight = rightValue;

    twoFragment->fLeft = oneFragment->fLeft;
    zeroFragment->fRight = oneFragment->fRight;
    twoFragment->leftPointTreeRefs[divDim] = oneFragment->leftPointTreeRefs[divDim];
    zeroFragment->rightPointTreeRefs[divDim] = oneFragment->rightPointTreeRefs[divDim];
}


template <std::size_t N>
void NParallel::Plane<N>::updateLambdaMax(double f1, double f2, double diff) {
    double newLambda = std::abs(f1 - f2) / diff;
    if (lambdaMax < newLambda) {
    #pragma omp critical
        if (lambdaMax < newLambda) {
            lambdaMax = newLambda;
        }
    }
}

template<std::size_t N>
void NParallel::Plane<N>::DivideFragment(
    NParallel::Plane<N>::NodePtr& fragment
) {
    searchFragments.InsertAfter(fragment, std::move(NParallel::Fragment<N>(*fragment->value)));
    fragment->next->value->Divide('0');
    searchFragments.InsertAfter(fragment->next, std::move(NParallel::Fragment<N>(*fragment->value)));
    fragment->next->next->value->Divide('2');
    fragment->value->Divide('1');
    
    tryAddFragmentPoints(fragment);
    
    double f1 = fragment->next->value->fRight;
    double f2 = fragment->next->value->fLeft;
    updateLambdaMax(f1, f2, fragment->value->diff);

    f1 = fragment->next->next->value->fRight;
    f2 = fragment->next->next->value->fLeft;
    updateLambdaMax(f1, f2, fragment->value->diff);

    f1 = fragment->value->fRight;
    f2 = fragment->value->fLeft;
    updateLambdaMax(f1, f2, fragment->value->diff);

    double diff = 2. * fragment->value->prevDiff / 3.;
    f1 = fragment->next->value->fLeft;
    f2 = fragment->value->fLeft;
    updateLambdaMax(f1, f2, diff);

    f1 = fragment->next->next->value->fRight;
    f2 = fragment->value->fRight;
    updateLambdaMax(f1, f2, diff);

    RecalcAfterDivision(fragment);
}


template <std::size_t N>
double NParallel::Plane<N>::GetBestPoint(std::array<double, N>& optimum) {
    return codeToPointId.GetBestPoint(optimum);
}


template <std::size_t N>
bool NParallel::Plane<N>::GetBestFragment(
    NParallel::Plane<N>::NodePtr& to
) {
    std::pair<double, NParallel::Plane<N>::NodePtr> topVal;
    bool resp = fragmentQueue.try_pop(topVal);
    if (resp) {
        to = topVal.second;
    }
    return resp;
}


template <std::size_t N>
void NParallel::Plane<N>::PutFragmentBack(
    NParallel::Plane<N>::NodePtr& fragment
) {
    // std::size_t k = std::max(searchFragments.size() / 2, 1ul);
    // double mu = (r + C / k) * lambdaMax;
    // double fLeft = fragment->value->fLeft;
    // double fRight = fragment->value->fRight;
    // fragment->value->UpdR(mu, fLeft, fRight);
    fragmentQueue.push(std::make_pair(fragment->value->R, fragment));
}


template <std::size_t N>
void NParallel::Plane<N>::RecalcAfterDivision(NParallel::Plane<N>::NodePtr& fragment) {
    tbb::concurrent_priority_queue<std::pair<double, NParallel::Plane<N>::NodePtr>> newFragmentQueue;
    std::size_t k = std::max(searchFragments.Size() / 2, 1ul);
    double mu = (r + C / k) * lambdaMax;

    std::size_t searchFragmentsSize = 3;
    while (searchFragmentsSize--) {
        double fLeft = fragment->value->fLeft;
        double fRight = fragment->value->fRight;
        fragment->value->UpdR(mu, fLeft, fRight);
        fragmentQueue.push(std::make_pair(fragment->value->R, fragment));
        fragment = fragment->next;
    }
}


template <std::size_t N>
void NParallel::Plane<N>::RecalcAllFragments() {
    tbb::concurrent_priority_queue<std::pair<double, NParallel::Plane<N>::NodePtr>> newFragmentQueue;
    std::size_t k = std::max(searchFragments.Size() / 2, 1ul);
    double mu = (r + C / k) * lambdaMax;
    prevMu = mu;

    std::size_t searchFragmentsSize = searchFragments.Size();
    NParallel::Plane<N>::NodePtr cur;
    searchFragments.GetBegin(cur);
    // bestFragment = cur;
    while (searchFragmentsSize--) {
        // if (AssertFragment(*cur->value)) {
        //     exit(0);
        // }
        double fLeft = cur->value->fLeft;
        double fRight = cur->value->fRight;
        cur->value->UpdR(mu, fLeft, fRight);
        // if (bestFragment->value->R < cur->value->R) {
        //     bestFragment = cur;
        // }
        newFragmentQueue.push(std::make_pair(cur->value->R, cur));
        cur = cur->next;
    }
    fragmentQueue.swap(newFragmentQueue);
}


