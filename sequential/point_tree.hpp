#include "one_dim_point_tree.hpp"

#include <array>
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>


namespace NSequential {

template <std::size_t N>
class PointTree {
public:
    std::array<double, N> left;
    std::array<double, N> right;
    std::function<double(const std::array<double, N>&)> f;

private:
    std::array<OneDimPointTree, N> oneDims;
    std::vector<double> pointValues;
    std::size_t bestPointId;
    std::array<double, N> bestPoint;

public:
    PointTree(
        const std::array<double, N>& left_, const std::array<double, N>& right_,
        const std::function<double(const std::array<double, N>&)>& f_)
        : left(left_)
        , right(right_)
        , f(f_) {
        for (OneDimPointTree& tree : oneDims) {
            tree.InitAllEdges();
        }
        bestPointId = 0;
        std::array<std::string, N> code;
        code.fill("0"); 
        AddPoint(code);
        code.fill("1");
        AddPoint(code);
    }

    std::size_t GetPointByCode(const std::array<std::string, N>& pointCode) {
        // std::cout << "\nGetPointByCode\n";
        std::map<std::size_t, std::size_t> pointIds;
        std::size_t candidate = 0;
        for (std::size_t i = 0; i != N; ++i) {
            // std::cout << pointCode[i] << std::endl;
            std::shared_ptr<OneDimPointTree> cur;
            oneDims[i].FindPointByCode(pointCode[i], cur);
            candidate = 0;
            for (const std::size_t& pointId : cur->pointIds) {
                // std::cout << pointId << ' ';
                if (i == pointIds[pointId]) {
                    candidate = pointId;
                    pointIds[pointId] = i + 1;
                }
            }
            // std::cout << std::endl;
            if (candidate == 0) {
                // std::cout << std::endl;
                return candidate;
            }
        }
        // std::cout << std::endl;
        return candidate;
    }

    std::size_t AddPoint(const std::array<std::string, N>& pointCode) {
        std::size_t pointId = pointValues.size()+1;
        std::array<double, N> x;
        // std::cout << pointId << ' ' << "f(";
        for (std::size_t i = 0; i != N; ++i) {
            // std::cout << '[' << pointCode[i] << "] ";
            double w = (pointCode[i][0] - '0');
            x[i] = left[i] + w * (right[i] - left[i]);
            // std::cout << "[" << x[i] << ", ";
            double p = 1.;
            for (char c : pointCode[i]) {
                if (p < 2) {
                    p *= 3.;
                    continue;
                }
                w = (c - '0');
                x[i] += w * (right[i] - left[i]) / p;
                // std::cout << x[i] << ", ";
                p *= 3.;
            }
            // std::cout << "]";
            // std::cout << x[i] << ' ';
            oneDims[i].AddPointCode(pointCode[i], pointId);
        }
        // std::cout << std::endl;
        pointValues.push_back(f(x));
        // std::cout << ") = " << pointValues.back() << std::endl;
        if (pointValues[bestPointId] > pointValues.back()) {
            bestPointId = pointId;
            std::swap(bestPoint, x);
        }
        return pointId;
    }

    void AddPointCode(const std::array<std::string, N>& pointCode, std::size_t pointId) {
        for (std::size_t i = 0; i != N; ++i) {
            oneDims[i].AddPointCode(pointCode[i], pointId);
        }
    }

    double GetValue(std::size_t pointId) const {
        assert(pointId > 0);
        return pointValues[pointId - 1];
    }

    std::size_t FCount() const {
        return pointValues.size();
    }

    double GetBestPoint(std::array<double, N>& optimum) const {
        optimum = bestPoint;
        return pointValues[bestPointId - 1];
    }
};
}