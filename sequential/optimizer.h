#pragma once
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>


namespace NSequential {

/**
 * Diagonal sequential optimizer for minimization
 * @tparam N dimensionality of optimized function
 */
template<std::size_t N>
class Optimizer {
private:

public:
    Optimizer();

    /**
     * Optimization will be done on this space:
     * @f$X = \{x | \forall i = 1 \cdots n \, left[i] \leq x[i] \leq right[i]\}@f$
     * @param f Lipschitz function to optimize
     * @param left Left point of hypercube optimization space
     * @param right Right point of hypercube optimization space
     */
    double optimize(
        const std::array<double, N>& left,
        const std::array<double, N>& right,
        const std::function<double(const std::array<double, N>&)>& f,
        std::array<double, N>& optimum,
        double eps = 0.02,
        double r = 2.8,
        double C = 100
    );
};

}


#include "optimizer.tpp"