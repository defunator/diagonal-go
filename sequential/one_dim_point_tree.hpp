#include <array>
#include <memory>
#include <string>
#include <cassert>
#include <vector>


namespace NSequential {

template <std::size_t N>
class PointTree;

class OneDimPointTree {
protected:
    std::vector<std::size_t> pointIds;

private:
    std::array<std::shared_ptr<OneDimPointTree>, 3> edges;

    std::size_t getLastNonZeroInd(const std::string& pointCode) {
        std::size_t res = pointCode.size() - 1;
        for (; res != 0 && pointCode[res] == '0'; --res);
        return res + 1;
    }

public:
    OneDimPointTree() { }

    template <std::size_t N>
    friend class PointTree;

    void InitAllEdges() {
        for (std::size_t i = 0; i != 3; ++i) {
            edges[i] = std::make_shared<OneDimPointTree>();
        }
    }

    void InitOneEdge(std::size_t id) {
        edges[id] = std::make_shared<OneDimPointTree>();
    }

    void AddPointId(std::size_t pointId) {
        pointIds.push_back(pointId);
    }

    void AddPointCode(const std::string& pointCode, std::size_t pointId) {
        // std::cout << "\nAddPointCode" << pointId << "\n";
        assert(pointCode.size() != 0);
        std::shared_ptr<OneDimPointTree> cur(edges[pointCode[0] - '0']);
        std::size_t lastNonZeroInd = getLastNonZeroInd(pointCode);
        // std::cout << lastNonZeroInd << ' ' << pointCode[0];

        for (std::size_t i = 1; i != lastNonZeroInd; ++i) {
            std::size_t id = pointCode[i] - '0';
            // std::cout << id << ' ';
            if (cur->edges[id] == nullptr) {
                cur->InitOneEdge(id);
            }
            cur = cur->edges[id];
        }
        // std::cout << "Added\n\n";
        cur->AddPointId(pointId);
    }

    void FindPointByCode(const std::string& pointCode, std::shared_ptr<OneDimPointTree>& result) {
        // std::cout << "\nFindPointByCode" << "\n";
        assert(pointCode.size() != 0);
        std::shared_ptr<OneDimPointTree> cur(edges[pointCode[0] - '0']);
        std::size_t lastNonZeroInd = getLastNonZeroInd(pointCode);
        // std::cout << lastNonZeroInd << ' ' << pointCode[0];

        for (std::size_t i = 1; i != lastNonZeroInd; ++i) {
            std::size_t id = pointCode[i] - '0';
            if (cur->edges[id] == nullptr) {
                cur->InitOneEdge(id);
            }
            cur = cur->edges[id];
            // std::cout << ' ' << id;
        }
        result = cur;
        // std::cout << "\n";
        // for (auto el : result->edges[0]->pointIds) {
        //     std::cout << el << ' ';
        // }
        // std::cout << "FindPointByCode\n\n";
    }
};

} // end NSequential namespace
