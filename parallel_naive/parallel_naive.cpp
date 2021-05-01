#include "optimizer.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>


using namespace std;


random_device rand_dev;
mt19937 generator(rand_dev());
uniform_real_distribution<double> distribution(0, 1.0);
const size_t N = 6;
const double relEps = 0.1;
// double fOptim = -0.352386; 
// double fOptim = -0.516374; 
// double fOptim = 1;
// double fOptim = 1;
double fOptim = 1;


double f(const array<double, N>& x) {
    // double x1 = x[0];
    // double x2 = x[1];
    // double x3 = x[2];

    // return 0.25 * x1 * x1 * x1 * x1 - 0.5 * x1 * x1 + 0.1 * x1 + 0.5 * x2 * x2; // [-10, 10]; min ~ -0.352386
    // return (x1 * x1 - 2 * x2 * x2 + x3 * x3) * sin(x1) * sin(x2) * sin(x3); // [-1, 1]; min ~ -0.516374
    
    // {
    //     array<double, N> y;
    //     for (size_t i = 0; i != N; ++i) {
    //         y[i] = 1 + (x[i] - 1) / 4;
    //     }
    //     double f = 0;
    //     f += 10 * sin(M_PI * y[0]) * sin(M_PI * y[0]);
    //     f += (y[N-1] - 1) * (y[N-1] - 1);
    //     for (size_t i = 0; i != N-1; ++i) {
    //         f += (y[i] - 1) * (y[i] - 1) * (1 + 10 * sin(M_PI * y[i + 1]) * sin(M_PI * y[i + 1]));
    //     }
    //     return M_PI * f / 10 + 1; // [-10, 10]; min = 1
    // }

    // return 1 + 100 * (x3 - 0.25 * (x1 + x2) * (x1 + x2)) * (x3 - 0.25 * (x1 + x2) * (x1 + x2)) + (1 - x1) * (1 - x1) + (1 - x2) * (1 - x2); // [0, 2] min = 0

    {
        double f = 0;
        f += sin(3 * M_PI * x[0]) * sin(3 * M_PI * x[0]);
        for (size_t i = 1; i != N - 1; ++i) {
            f += (x[i] - 1) * (x[i] - 1) * (1 + sin(3 * M_PI * x[i + 1]) * sin(3 * M_PI * x[i + 1]));
        }
        f += (x[N-1] - 1) * (x[N-1] - 1) * (1 + sin(2 * M_PI * x[N-1]) * sin(2 * M_PI * x[N-1]));
        f *= 0.1;
        return f + 1; // [-5, 5]; min = 1
    }
}

double getScore(double fSearch) {
    return std::abs(fSearch - fOptim) / std::abs(fOptim);
}


double monteCarloSearch(const array<double, N>& left,
                const array<double, N>& right, array<double, N>& res) {
    res = left;
    size_t f_count = 0;
    while (42) {
        array<double, N> pt;
        for (size_t i = 0; i != N; ++i) {
            pt[i] = left[i] + distribution(generator) * (right[i] - left[i]);
        }
        if (f(pt) < f(res)) {
            res = std::move(pt);
        }
        if (getScore(f(res)) < relEps) {
            break;
        }
        ++f_count;
        if (f_count % 1000000 == 0) {
            cout << f_count << ' ' << getScore(f(res)) << endl;
        }
    }
    cout << "F count = " << f_count << endl;
    return f(res);
}


int main() {
    NSequential::Optimizer<N> op;
    array<double, N> left;
    array<double, N> right;
    left.fill(-5);
    right.fill(5);
    array<double, N> res;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    double optimum = op.optimize(left, right, f, res, fOptim);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "DIAGONAL RESULT:" << endl;
    cout << "X = (";
    for (size_t i = 0; i != res.size(); ++i) {
        cout << res[i];
        if (i != res.size() - 1) {
            cout << ", ";
        }
    }
    cout << ")" << endl;
    cout << "f(X) = " << optimum << endl;
    cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl;

    // cout << "BRUTE FORCE RESULT:" << endl;
    // begin = chrono::steady_clock::now();
    // double realOpt = bruteSearch(left, right, res);
    // end = chrono::steady_clock::now();
    // cout << "X = (";
    // for (size_t i = 0; i != res.size(); ++i) {
    //     cout << res[i];
    //     if (i != res.size() - 1) {
    //         cout << ", ";
    //     }
    // }
    // cout << ")" << endl;
    // cout << "f(X) = " << realOpt << endl;
    // cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl;

    cout << "MONTE-CARLO RESULT:" << endl;
    begin = chrono::steady_clock::now();
    double monteCarlo = monteCarloSearch(left, right, res);
    end = chrono::steady_clock::now();
    cout << "X = (";
    for (size_t i = 0; i != res.size(); ++i) {
        cout << res[i];
        if (i != res.size() - 1) {
            cout << ", ";
        }
    }
    cout << ")" << endl;
    cout << "f(X) = " << monteCarlo << endl;
    cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl;

    return 0;
}
