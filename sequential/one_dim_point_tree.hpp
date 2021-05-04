#include <array>
#include <memory>
#include <string>
#include <cassert>
#include <vector>


namespace NSequential {

template <std::size_t N>
class PointTree;

class OneDimPointTree {
public:
    double x;
    double d;
private:
    std::vector<std::size_t> pointIds;

private:
    std::array<std::shared_ptr<OneDimPointTree>, 3> edges;

    std::size_t getLastNonZeroInd(const std::string& pointCode) {
        std::size_t res = pointCode.size() - 1;
        for (; res != 0 && pointCode[res] == '0'; --res);
        return res + 1;
    }

public:
    OneDimPointTree(double x_, double d_)
        : x(x_)
        , d(d_) { }

    template <std::size_t N>
    friend class PointTree;

    void InitAllEdges() {
        for (std::size_t i = 0; i != 3; ++i) {
            edges[i] = std::make_shared<OneDimPointTree>(x + d * i, d / 3.);
        }
    }

    void InitOneEdge(std::size_t id) {
        edges[id] = std::make_shared<OneDimPointTree>(x + d * id, d / 3.);
    }

    void AddPointId(std::size_t pointId) {
        pointIds.push_back(pointId);
    }

    void AddPointCode(const std::string& pointCode, std::size_t pointId) {
        if (pointCode.empty()) {
            AddPointId(pointId);
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
        cur->AddPointId(pointId);
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
};

} // end NSequential namespace
