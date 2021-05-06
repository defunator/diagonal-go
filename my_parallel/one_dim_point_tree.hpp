#pragma once
#include "atomic_list.hpp"

#include <array>
#include <atomic>
#include <memory>
#include <string>
#include <cassert>
#include <vector>


namespace NParallel {

template <std::size_t N>
class PointTree;

class OneDimPointTree {
public:
    double x;
    double d;
private:
    AtomicList<std::pair<std::size_t, double>> points;

    enum EdgeStatus {
        EMPTY,
        WORKING,
        FILLED
    };

private:
    std::array<std::shared_ptr<OneDimPointTree>, 3> edges;
    std::array<std::atomic<EdgeStatus>, 3> edgeStatus;

    std::size_t getLastNonZeroInd(const std::string& pointCode) {
        std::size_t res = pointCode.size() - 1;
        for (; res != 0 && pointCode[res] == '0'; --res);
        return res + 1;
    }

public:
    OneDimPointTree(double x_, double d_)
        : x(x_)
        , d(d_) {
        for (auto& el : edgeStatus) {
            el = EdgeStatus::EMPTY;
        }
    }

    template <std::size_t N>
    friend class PointTree;

    void InitOneEdge(std::size_t id) {
        EdgeStatus curStatus{EdgeStatus::EMPTY};
        if (edgeStatus[id].compare_exchange_strong(curStatus, EdgeStatus::WORKING)) {
            edges[id] = std::make_shared<OneDimPointTree>(x + d * id, d / 3.);
            curStatus = EdgeStatus::WORKING;
            assert(edgeStatus[id].compare_exchange_strong(curStatus, EdgeStatus::FILLED));
        } else if (curStatus == EdgeStatus::WORKING) {
            curStatus = EdgeStatus::FILLED;
            while (!edgeStatus[id].compare_exchange_weak(curStatus, EdgeStatus::FILLED)) {
                curStatus = EdgeStatus::FILLED;
            }
        }
    }

    void InitAllEdges() {
        for (std::size_t i = 0; i != 3; ++i) {
            InitOneEdge(i);
        }
    }

    void AddPoint(std::size_t pointId, double pointValue) {
        points.InsertBack(std::make_pair(pointId, pointValue));
    }

    void AddPointCode(const std::string& pointCode, std::size_t pointId, double pointValue) {
        if (pointCode.empty()) {
            AddPoint(pointId, pointValue);
            return;
        }
        std::shared_ptr<OneDimPointTree> cur(edges[pointCode[0] - '0']);
        std::size_t lastNonZeroInd = getLastNonZeroInd(pointCode);

        for (std::size_t i = 1; i != lastNonZeroInd; ++i) {
            std::size_t id = pointCode[i] - '0';
            if (cur->edges[id] == nullptr) {
                cur->InitOneEdge(id);
            }
            cur = cur->edges[id];
        }
        cur->AddPoint(pointId, pointValue);
    }

    void FindPointByCode(const std::string& pointCode, std::shared_ptr<OneDimPointTree>& result) {
        std::shared_ptr<OneDimPointTree> cur(edges[pointCode[0] - '0']);
        std::size_t lastNonZeroInd = getLastNonZeroInd(pointCode);

        for (std::size_t i = 1; i != lastNonZeroInd; ++i) {
            std::size_t id = pointCode[i] - '0';
            if (cur->edges[id] == nullptr) {
                cur->InitOneEdge(id);
            }
            cur = cur->edges[id];
        }
        result = cur;
    }

    void Traverse(char planeId, std::shared_ptr<OneDimPointTree>& to) {
        std::size_t id = planeId - '0';
        if (edges[id] == nullptr) {
            InitOneEdge(id);
        }
        to = edges[id];
    }

    void Traverse(char planeId, std::size_t zeroDivsCount, std::shared_ptr<OneDimPointTree>& to) {
        bool first = true;
        std::shared_ptr<OneDimPointTree> tmp;
        while (zeroDivsCount--) {
            if (first) {
                Traverse('0', tmp);
                first = false;
            } else {
                to->Traverse('0', tmp);
            }
            to = tmp;
        }
        if (first) {
            Traverse(planeId, tmp);
            first = false;
        } else {
            to->Traverse(planeId, tmp);
        }
        to = tmp;
    }

    void dfs() {
        for (auto& el : edges) {
            if (el) {
                el->dfs();
                el.reset();
            }
        }
    }

    ~OneDimPointTree() {
        dfs();
    }
};

} // end NParallel namespace
