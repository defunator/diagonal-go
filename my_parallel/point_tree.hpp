#pragma once
#include "atomic_list.hpp"
#include "one_dim_point_tree.hpp"

#include <array>
#include <atomic>
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <omp.h>
#include <queue>
#include <string>
#include <tbb/concurrent_vector.h>
#include <unordered_map>
#include <vector>


namespace NParallel {

template <std::size_t N>
class PointTree {
public:
    std::array<double, N> left;
    std::array<double, N> right;
    std::function<double(const std::array<double, N>&)> f;

public:
    std::atomic<std::size_t> fCount;
    std::array<std::shared_ptr<OneDimPointTree>, N> oneDims;
    // AtomicList<double> pointValues;
    std::atomic<double> bestValue;
    std::atomic<std::array<double, N>> bestValuePoint;

public:
    PointTree(
        const std::array<double, N>& left_, const std::array<double, N>& right_,
        const std::function<double(const std::array<double, N>&)>& f_)
        : left(left_)
        , right(right_)
        , f(f_)
        , fCount(0) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i] = std::make_shared<OneDimPointTree>(left[i], right[i] - left[i]);
            oneDims[i]->InitAllEdges();
        }
        std::array<std::string, N> code;
        code.fill("0"); 
        AddPoint(code);
        code.fill("1");
        AddPoint(code);
    }

    // std::size_t GetPointByCode(const std::array<std::string, N>& pointCode) {
    //     std::map<std::size_t, std::size_t> pointIds;
    //     std::size_t candidate = 0;
    //     for (std::size_t i = 0; i != N; ++i) {
    //         std::shared_ptr<OneDimPointTree> cur;
    //         oneDims[i]->FindPointByCode(pointCode[i], cur);
    //         candidate = 0;
    //         for (const std::size_t& pointId : cur->pointIds) {
    //             if (i == pointIds[pointId]) {
    //                 candidate = pointId;
    //                 pointIds[pointId] = i + 1;
    //             }
    //         }
    //         if (candidate == 0) {
    //             return candidate;
    //         }
    //     }
    //     return candidate;
    // }

    std::pair<std::size_t, double> GetPointIdFromArray(const std::array<std::shared_ptr<OneDimPointTree>, N>& pointNodes) const {
        std::priority_queue<std::pair<std::size_t, std::size_t>> pointNodesQueue;
        for (std::size_t i = 0; i != N; ++i) {
            if (pointNodes[i]->points.Size() == 0) {
                return {0, 0};
            }
            pointNodesQueue.emplace(-pointNodes[i]->points.Size(), i);
        }

        std::unordered_map<std::size_t, std::size_t> pointIds;
        std::size_t candidateId = 0;
        double candidateValue;
        std::size_t it_num = 0;
        while (!pointNodesQueue.empty()) {
            std::size_t i = pointNodesQueue.top().second;
            candidateId = 0;

            std::size_t listSize = pointNodes[i]->points.Size();
            std::shared_ptr<Node<std::pair<std::size_t, double>>> it;
            pointNodes[i]->points.GetBegin(it);

            while (listSize--) {
                std::size_t pointId = it->value->first;
                auto tmp = pointIds.find(pointId);
                if (it_num == 0 || (tmp != pointIds.end() && it_num == tmp->second)) {
                    candidateId = pointId;
                    candidateValue = it->value->second;
                    pointIds[pointId] = it_num + 1;
                }
                it = it->next;
            }

            if (candidateId == 0) {
                return {candidateId, 0};
            }
            pointNodesQueue.pop();
            ++it_num;
        }
        return {candidateId, candidateValue};
    }

    void GetTreeNodesByCode(const std::array<std::string, N>& pointCode,
                            std::array<std::shared_ptr<OneDimPointTree>, N>& pointNodes) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i]->FindPointByCode(pointCode[i], pointNodes[i]);
        }
    }

    // void GetTreeNodeByCode(const std::string& pointCode,
    //                         std::shared_ptr<OneDimPointTree>& pointNode,
    //                         std::size_t dim) {
    //     oneDims[dim]->FindPointByCode(pointCode, pointNode);
    // }

    std::pair<std::size_t, double> AddPoint(const std::array<std::string, N>& pointCode) {
        std::size_t pointId = fCount.fetch_add(1) + 1;
        std::array<double, N> x;
        for (std::size_t i = 0; i != N; ++i) {
            double w = (pointCode[i][0] - '0');
            x[i] = left[i] + w * (right[i] - left[i]);
            double p = 1.;
            for (char c : pointCode[i]) {
                if (p < 2) {
                    p *= 3.;
                    continue;
                }
                w = (c - '0');
                x[i] += w * (right[i] - left[i]) / p;
                p *= 3.;
            }
        }
        double pointValue = f(x);
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i]->AddPointCode(pointCode[i], pointId, pointValue);
        }
        if (pointId == 1 || bestValue > pointValue) {
        #pragma omp critical
            if (pointId == 1 || bestValue > pointValue) {
                bestValue = pointValue;
                bestValuePoint = std::move(x);
            }
        }
        return {pointId, pointValue};
    }

    std::pair<std::size_t, double> AddPoint(std::array<std::shared_ptr<OneDimPointTree>, N>& pointTreeRefs) {
        std::size_t pointId = fCount.fetch_add(1) + 1;
        std::array<double, N> x;
        for (std::size_t i = 0; i != N; ++i) {
            x[i] = pointTreeRefs[i]->x;
        }
        double pointValue = f(x);
        for (std::size_t i = 0; i != N; ++i) {
            pointTreeRefs[i]->AddPointCode("", pointId, pointValue);
        }

        if (bestValue > pointValue) {
        #pragma omp critical
            if (bestValue > pointValue) {
                bestValue = pointValue;
                bestValuePoint = std::move(x);
            }
        }
        return {pointId, pointValue};
    }

    // void AddPointCode(const std::array<std::string, N>& pointCode, std::size_t pointId) {
    //     for (std::size_t i = 0; i != N; ++i) {
    //         oneDims[i]->AddPointCode(pointCode[i], pointId);
    //     }
    // }

    // void AddPointCode(std::array<std::shared_ptr<OneDimPointTree>, N>& pointTreeRefs, std::size_t pointId) {
    //     for (std::size_t i = 0; i != N; ++i) {
    //         pointTreeRefs[i]->AddPointCode("", pointId);
    //     }
    // }

    // double GetValue(std::size_t pointId) const {
    //     Node<double> tmp;
    //     pointValues.GetBegin(tmp);
    //     while (--pointId) {
    //         tmp = tmp->next;
    //     }
    //     return tmp->value;
    // }

    std::size_t FCount() const {
        return fCount;
    }

    double GetBestPoint(std::array<double, N>& optimum) const {
        optimum = bestValuePoint;
        return bestValue;
    }

    double GetBestPoint() const {
        return bestValue;
    }

    // ~PointTree() {
    //     for (std::size_t i = 0; i != N; ++i) {
    //         oneDims[i]->dfs();
    //     }
    // }
};
}