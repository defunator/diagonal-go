#include "one_dim_point_tree.hpp"

#include <array>
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>


namespace NSequential {

template <std::size_t N>
class PointTree {
public:
    std::array<double, N> left;
    std::array<double, N> right;
    std::function<double(const std::array<double, N>&)> f;

public:
    std::array<std::shared_ptr<OneDimPointTree>, N> oneDims;
    std::vector<double> pointValues;
    // std::vector<std::array<double, N>> points;
    std::size_t bestPointId;
    std::array<double, N> bestPoint;

public:
    PointTree(
        const std::array<double, N>& left_, const std::array<double, N>& right_,
        const std::function<double(const std::array<double, N>&)>& f_)
        : left(left_)
        , right(right_)
        , f(f_) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i] = std::make_shared<OneDimPointTree>(left[i], right[i] - left[i]);
            oneDims[i]->InitAllEdges();
        }
        pointValues.reserve(20000);
        // points.reserve(20000);
        bestPointId = 0;
        std::array<std::string, N> code;
        code.fill("0"); 
        AddPoint(code);
        code.fill("1");
        AddPoint(code);
    }

    std::size_t GetPointByCode(const std::array<std::string, N>& pointCode) {
        std::map<std::size_t, std::size_t> pointIds;
        std::size_t candidate = 0;
        for (std::size_t i = 0; i != N; ++i) {
            std::shared_ptr<OneDimPointTree> cur;
            oneDims[i]->FindPointByCode(pointCode[i], cur);
            candidate = 0;
            for (const std::size_t& pointId : cur->pointIds) {
                if (i == pointIds[pointId]) {
                    candidate = pointId;
                    pointIds[pointId] = i + 1;
                }
            }
            if (candidate == 0) {
                return candidate;
            }
        }
        return candidate;
    }

    std::size_t GetPointIdFromArray(const std::array<std::shared_ptr<OneDimPointTree>, N>& pointNodes) const {
        std::priority_queue<std::pair<std::size_t, std::size_t>> pointNodesQueue;
        for (std::size_t i = 0; i != N; ++i) {
            if (pointNodes[i]->pointIds.empty()) {
                return 0;
            }
            pointNodesQueue.emplace(-pointNodes[i]->pointIds.size(), i);
        }

        std::unordered_map<std::size_t, std::size_t> pointIds;
        std::size_t candidate = 0;
        std::size_t it_num = 0;
        while (!pointNodesQueue.empty()) {
            std::size_t i = pointNodesQueue.top().second;
            candidate = 0;
            for (const std::size_t& pointId : pointNodes[i]->pointIds) {
                auto it = pointIds.find(pointId);
                if (it_num == 0 || (it != pointIds.end() && it_num == it->second)) {
                    candidate = pointId;
                    pointIds[pointId] = it_num + 1;
                }
            }
            if (candidate == 0) {
                return candidate;
            }
            pointNodesQueue.pop();
            ++it_num;
        }
        return candidate;
    }

    void GetTreeNodesByCode(const std::array<std::string, N>& pointCode,
                            std::array<std::shared_ptr<OneDimPointTree>, N>& pointNodes) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i]->FindPointByCode(pointCode[i], pointNodes[i]);
        }
    }

    void GetTreeNodeByCode(const std::string& pointCode,
                            std::shared_ptr<OneDimPointTree>& pointNode,
                            std::size_t dim) {
        oneDims[dim]->FindPointByCode(pointCode, pointNode);
    }

    std::size_t AddPoint(const std::array<std::string, N>& pointCode) {
        std::size_t pointId = pointValues.size()+1;
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
            oneDims[i]->AddPointCode(pointCode[i], pointId);
        }
        pointValues.push_back(f(x));
        // points.push_back(std::move(x));
        if (pointValues[bestPointId] > pointValues.back()) {
            bestPointId = pointId-1;
            // bestPoint = points.back();
            bestPoint = std::move(x);
        }
        return pointId;
    }

    // std::size_t AddPoint(std::size_t leftPointId, std::size_t rightPointId, std::size_t dim, bool leftP = true) {
    //     std::array<double, N> x(points[leftPointId - 1]);
    //     if (!leftP) {
    //         x = points[rightPointId - 1];
    //         x[dim] += 2. * (points[leftPointId - 1][dim] - points[rightPointId - 1][dim]) / 3.;
    //     } else {
    //         x[dim] += 2. * (points[rightPointId - 1][dim] - points[leftPointId - 1][dim]) / 3.;
    //     }
    //     std::size_t pointId = pointValues.size() + 1;
    //     pointValues.push_back(f(x));
    //     // points.push_back(std::move(x));
    //     if (pointValues[bestPointId] > pointValues.back()) {
    //         bestPointId = pointId-1;
    //         // bestPoint = points.back();
    //         bestPoint = std::move(x);
    //     }
    //     return pointId;
    // }

    std::size_t AddPoint(std::array<std::shared_ptr<OneDimPointTree>, N>& pointTreeRefs) {
        std::size_t pointId = pointValues.size() + 1;
        std::array<double, N> x;
        for (std::size_t i = 0; i != N; ++i) {
            x[i] = pointTreeRefs[i]->x;
            pointTreeRefs[i]->AddPointCode("", pointId);
        }
        pointValues.push_back(f(x));
        // points.push_back(std::move(x));
        if (pointValues[bestPointId] > pointValues.back()) {
            bestPointId = pointId-1;
            // bestPoint = points.back();
            bestPoint = std::move(x);
        }
        return pointId;
    }

    void AddPointCode(const std::array<std::string, N>& pointCode, std::size_t pointId) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i]->AddPointCode(pointCode[i], pointId);
        }
    }

    void AddPointCode(std::array<std::shared_ptr<OneDimPointTree>, N>& pointTreeRefs, std::size_t pointId) {
        for (std::size_t i = 0; i != N; ++i) {
            pointTreeRefs[i]->AddPointCode("", pointId);
        }
    }

    double GetValue(std::size_t pointId) const {
        assert(pointId > 0);
        return pointValues[pointId - 1];
    }

    // double GetValue(std::size_t pointId, std::array<double, N>& x) const {
    //     assert(pointId > 0);
    //     x = points[pointId - 1];
    //     return pointValues[pointId - 1];
    // }

    std::size_t FCount() const {
        return pointValues.size();
    }

    double GetBestPoint(std::array<double, N>& optimum) const {
        optimum = bestPoint;
        return pointValues[bestPointId];
    }
};
}