// #include "optimizer.hpp"
// #include "parallel_naive.hpp"
#include "optimizer_concurrent_queue.hpp"
#include "tests.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>


using namespace std;

// const double relEps = 0.1;

template <std::size_t N>
void benchmarkTest(const NTest::Test<N>& test) {
    // NSequential::Optimizer<N> optimizer;
    // NParallel::Optimizer<N> optimizer;
    NParallel::Optimizer<N> optimizer;
    array<double, N> xMin;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    double optimum = optimizer.optimize(test, xMin);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << test.GetName() << endl;
    cout << "OPTIMUM: f(";
    for (std::size_t i = 0; i != N; ++i) {
        cout << test.GetArgMin()[i];
        if (i != N - 1) {
            cout << ", ";
        }
    }
    cout << ") = " << test.GetMinValue() << endl;

    cout << "DIAGONAL RESULT:" << endl;
    cout << "f(";
    for (size_t i = 0; i != N; ++i) {
        cout << xMin[i];
        if (i != N - 1) {
            cout << ", ";
        }
    }
    cout << ") = " << optimum << endl;
    cout << "TIME = " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << endl << endl;
}

int main() {
    // NParallel::Optimizer<N> op;
    // benchmarkTest(NTest::GetTest0());
    // benchmarkTest(NTest::GetTest1());
    // benchmarkTest(NTest::GetTest2());
    // benchmarkTest(NTest::GetTestChichinadze());
    // // // // benchmarkTest(NTest::GetTest4<6>());
    // benchmarkTest(NTest::GetTest5());
    // benchmarkTest(NTest::GetTest6<5>());
    benchmarkTest(NTest::GetTestHartman());
    // benchmarkTest(NTest::GetTest3());
    // benchmarkTest(NTest::GetTestDolan());
    

    return 0;
}
