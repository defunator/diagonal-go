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


double f(const array<double, 2>& x) {
    double x1 = x[0];
    double x2 = x[1];

    // return 0.25 * x1 * x1 * x1 * x1 - 0.5 * x1 * x1 + 0.1 * x1 + 0.5 * x2 * x2; // [-10, 10]
    // return -4 * x1 * x2 * sin(4 * M_PI * x2); // [0, 1]
    return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (x1 - 1) * (x1 - 1); // [-2, 8]
    // return (x1 * x1 + x2 - 11) * (x1 * x1 + x2 - 11) + (x1 + x2 * x2 - 7) * (x1 + x2 * x2 - 7); // [-2, 6]
    // return 2 * x1 * x1 - 1.05 * x1 * x1 * x1 * x1 + x1 * x1 * x1 * x1 * x1 * x1 / 6 + x1 * x2 + x2 * x2; // [-5, 5]
}


double monteCarloSearch(const array<double, 2>& left,
                const array<double, 2>& right, array<double, 2>& res) {
    double d = 1e-4;
    res = left;
    size_t f_count = 0;
    while (42) {
        array<double, 2> tmp_res(res);
        for (size_t i = 0; i != 400; ++i, ++f_count) {
            array<double, 2> pt;
            pt[0] = left[0] + distribution(generator) * (right[0] - left[0]);
            pt[1] = left[1] + distribution(generator) * (right[1] - left[1]);
            if (f(pt) < f(tmp_res)) {
                tmp_res = pt;
            }
        }
        if (f(res) - f(tmp_res) < d) {
            break;
        }
        res = tmp_res;
    }
    cout << "F count = " << f_count << endl;
    return f(res);
}


int main() {
    NSequential::Optimizer<2> op;
    array<double, 2> left;
    array<double, 2> right;
    left.fill(-2);
    right.fill(8);
    array<double, 2> res;

    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    double optimum = op.optimize(left, right, f, res);
    // // chrono::steady_clock::time_point end = chrono::steady_clock::now();

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
    // cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl;


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
    // begin = chrono::steady_clock::now();
    double monteCarlo = monteCarloSearch(left, right, res);
    // end = chrono::steady_clock::now();
    cout << "X = (";
    for (size_t i = 0; i != res.size(); ++i) {
        cout << res[i];
        if (i != res.size() - 1) {
            cout << ", ";
        }
    }
    cout << ")" << endl;
    cout << "f(X) = " << monteCarlo << endl;
    // cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl;

    return 0;
}
