#include "optimizer.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <random>


using namespace std;


random_device rand_dev;
mt19937 generator(rand_dev());
uniform_real_distribution<double> distribution(-1.0, 1.0);
const size_t N_GKLS = 7;
vector<vector<double>> A(N_GKLS, vector<double>(N_GKLS));
vector<vector<double>> B(N_GKLS, vector<double>(N_GKLS));
vector<vector<double>> C(N_GKLS, vector<double>(N_GKLS));
vector<vector<double>> D(N_GKLS, vector<double>(N_GKLS));


void initGKLS() {
    for (size_t i = 0; i != N_GKLS; ++i) {
        for (size_t j = 0; j != N_GKLS; ++j) {
            A[i][j] = distribution(generator);
            B[i][j] = distribution(generator);
            C[i][j] = distribution(generator);
            D[i][j] = distribution(generator);
        }
    }
    // for (size_t i = 0; i != N_GKLS; ++i) {
    //     for (size_t j = 0; j != N_GKLS; ++j) {
    //         cout << A[i][j] << "* sin(pi*" << i << "*x_1)*sin(pi*" << j << "*x_2)+";
    //         cout << B[i][j] << "* cos(pi*" << i << "*x_1)*cos(pi*" << j << "*x_2)+";
    //     }
    // }
    // cout << endl;
    // for (size_t i = 0; i != N_GKLS; ++i) {
    //     for (size_t j = 0; j != N_GKLS; ++j) {
    //         cout << C[i][j] << "* sin(pi*" << i << "*x_1)*sin(pi*" << j << "*x_2)+";
    //         cout << D[i][j] << "* cos(pi*" << i << "*x_1)*cos(pi*" << j << "*x_2)+";
    //     }
    // }
    // cout << endl;
}


double f(const array<double, 2>& x) {
    // array<double, 10> y;
    // for (size_t i = 0; i != 10; ++i) {
    //     y[i] = 1 + (x[i] - 1) / 4;
    // }
    // double f = 0;
    // f += 10 * sin(M_PI * y[0]) * sin(M_PI * y[0]);
    // f += (y[9] - 1) * (y[9] - 1);
    // for (size_t i = 0; i != 9; ++i) {
    //     f += (y[i] - 1) * (y[i] - 1) * (1 + 10 * sin(M_PI * y[i + 1]) * sin(M_PI * y[i + 1]));
    // }
    // return M_PI * f / 10;

    double x1 = x[0];
    double x2 = x[1];
    // double x3 = x[2];

    // double s1 = 0;
    // double s2 = 0;
    // for (size_t i = 0; i != N_GKLS; ++i) {
    //     for (size_t j = 0; j != N_GKLS; ++j) {
    //         double a = sin((i + 1) * M_PI * x1) * sin((j + 1) * M_PI * x2);
    //         double b = cos((i + 1) * M_PI * x1) * cos((j + 1) * M_PI * x2);
    //         s1 += A[i][j] * a + B[i][j] * b;
    //         s2 += C[i][j] * a + D[i][j] * b;
    //     }
    // }
    // double f = -sqrt(s1 * s1 + s2 * s2);
    // return f;

    // return 0.01 * (x1 * x2 + (x1 - M_PI) * (x1 - M_PI) + 3 * (x2 - M_PI) * (x2 - M_PI)) - abs(sin(x1) * sin(2 * x2));
    // return (x1 * x1 - 2 * x2 * x2 + x3 * x3) * sin(x1) * sin(x2) * sin(x3);
    // return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (x1 - 1) * (x1 - 1);
    // return 2 * x1 * x1 - 1.05 * x1 * x1 * x1 * x1 + x1 * x1 * x1 * x1 * x1 * x1 / 6 + x1 * x2 + x2 * x2;
    // return -sin(2 * x1 + 1) - 2 * sin(3 * x2 + 2);
    // return -4 * x1 * x2 * sin(4 * M_PI * x2);
    return 0.25 * x1 * x1 * x1 * x1 - 0.5 * x1 * x1 + 0.1 * x1 + 0.5 * x2 * x2;
}

// double bruteSearch(const array<double, 2>& left, const array<double, 2>& right) {
//     double d = 1e-4;
//     array<double, 2> best({0, 0});
//     for (double x_1 = left[0]; x_1 < right[0]; x_1 += d) {
//         for (double x_2 = left[1]; x_2 < right[1]; x_2 += d) {
//             array<double, 2> x({x_1, x_2});
//             if (f(x) < f(best)) {
//                 best = x;
//             }
//         }
//     }
//     for (auto el : best) {
//         cout << el << ' ';
//     }
//     cout << endl;
//     return f(best);
// }


int main() {
    initGKLS();

    NSequential::Optimizer<2> op;
    array<double, 2> left;
    array<double, 2> right;
    left.fill(-10);
    right.fill(10);
    array<double, 2> res;
    double optimum = op.optimize(left, right, f, res);
    cout << "RESULT:" << endl;
    cout << "X = (";
    for (double el: res) {
        cout << el << ", ";
    }
    cout << ')' << endl;
    cout << "y = " << optimum << endl;
    // cout << "REAL BEST:" << endl;
    // double realOpt = bruteSearch(left, right);
    // cout << "f(X*) = " << realOpt << endl;
}
