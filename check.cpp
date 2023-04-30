#include<iostream>
#include<vector>
#include "matrix.h"
#include "rational.h"
#include "vector.h"

namespace Matrix_tests {
    void fail(const char *message) {
        std::cerr << "FAIL: ";
        std::cerr << message;
        exit(0);
    }
}
namespace Tasks {
    void Task1() {
        std::cout << "-----Task 1-----\n";
        Algebra::Matrix<double> A({{6,  3,  -1, 1},
                                   {-1, 1,  1,  -2},
                                   {-1, -4, 6,  -2},
                                   {1,  4,  -1, 7}});
        Algebra::Matrix<double> B({{4,  -7, -5, -3},
                                   {-1, 1,  -3, -2},
                                   {1,  4,  8,  2},
                                   {1,  4,  3,  7}});
        Algebra::Matrix<double> C({{6,  4,  2,  -1},
                                   {-1, -1, -3, 2},
                                   {1,  8,  9,  -3},
                                   {-1, -4, -2, 6}});
        std::cout << "A: " << A.Characteristic_poly() << "\n";
        std::cout << "B: " << B.Characteristic_poly() << "\n";
        std::cout << "C: " << C.Characteristic_poly() << "\n";
        std::cout << "rk(A - 5) = " << (A - 5).rank() << "\n";
        std::cout << "rk(B - 5) = " << (B - 5).rank() << "\n";
        std::cout << "rk(C - 5) = " << (C - 5).rank() << "\n";
        std::cout << "rk((A - 5)^2) = " << ((A - 5) * (A - 5)).rank() << "\n";
        std::cout << "rk((B - 5)^2) = " << ((B - 5) * (B - 5)).rank() << "\n";
        std::cout << (C - 5) * (C - 5) << "\n";
        std::cout << "rk((C - 5)^2) = " << ((C - 5) * (C - 5)).rank() << "\n";
    }

    void Task2() {
        std::cout << "-----Task 2-----\n";
        Algebra::Matrix<double> A({{3,  -3, -1, 0, 0,  0},
                                   {0,  5,  1,  0, 0,  0},
                                   {0,  -2, 2,  0, 0,  0},
                                   {-2, -3, -3, 1, -1, 0},
                                   {4,  8,  6,  4, 5,  0},
                                   {1,  2,  1,  2, 2,  4}});
        Algebra::Matrix<double> B({{3, -3, -1},
                                   {0, 5,  1},
                                   {0, -2, 2}});
        Algebra::Matrix<double> D({{1, -1, 0},
                                   {4, 5,  0},
                                   {2, 2,  4}});
        std::cout << "B : " << B.Characteristic_poly() << "\n";
        std::cout << "C : " << D.Characteristic_poly() << "\n";
        std::cout << A.Characteristic_poly() << "\n";
        std::cout << "rk(A - 3) = " << (A - 3).rank() << "\n";
        std::cout << "rk((A - 3)^2) = " << ((A - 3) * (A - 3)).rank() << "\n";
        std::cout << "rk((A - 3)^3) = " << ((A - 3) * (A - 3) * (A - 3)).rank() << "\n\n";
        std::cout << "rk(A - 4) = " << (A - 4).rank() << "\n";
        std::cout << "A = \n" << (A - 3).row_echelon_form() << "\n";
        Algebra::Matrix<double> u1({{1, 0, 0, -1, 0, 1}});
        Algebra::Matrix<double> u2({{-1, 0, 0, -1.5, 1, 0}});
        std::cout << "(A - 3)^2:\n" << ((A - 3) * (A - 3)).row_echelon_form() << "\n";
        Algebra::Matrix<double> v1({{-1, 0, 0, 0, 0, 5}});
        Algebra::Matrix<double> v2({{-4, 0, 0, 0, 5, 0}});
        Algebra::Matrix<double> v3({{-6, 0, 0, 5, 0, 0}});
        Algebra::Matrix<double> v4({{-3, -5, 10, 0, 0, 0}});
        std::cout << "f1 :\n" << ((A - 3) * v1.T()) << "\n";
        std::cout << "f3 :\n" << ((A - 3) * v4.T()) << "\n";
        std::cout << "(A - 4)^2 :\n" << ((A - 4) * (A - 4)).row_echelon_form() << "\n";
        Algebra::Matrix<double> r1({{0, 0, 0, 0, 0, 1}});
        Algebra::Matrix<double> r2({{2, -1, 1, -2, 2, 0}});
        std::cout << "f5 :\n" << (A - 4) * r2.T() << "\n";
        Algebra::Matrix<double> Base({{0,  -1, 0,  -3, 0, 2},
                                      {0,  0,  0,  -5, 0, -1},
                                      {0,  0,  0,  10, 0, 1},
                                      {2,  0,  -9, 0,  0, -2},
                                      {-4, 0,  8,  0,  0, 2},
                                      {4,  5,  -3, 0,  1, 0}});
        std::cout << Base.rank();
    }

    void Task3() {
        std::cout << "-----Task 3-----\n";
        Algebra::Matrix<double> g1({{-1, 2, 2}}), g2({{2, -3, -2}}), g3({{2, -2, -1}});
        Algebra::Matrix<double> c1({{-1, -2, 2}}), c2({{-2, -3, 2}}), c3({{2, 2, -1}});
        Algebra::Matrix<double> A({{-2, -3, -5},
                                   {-6, -3, -9},
                                   {6,  1,  7}});
        std::cout << A * c1.T() << "\n" << A * c2.T() << "\n" << A * c3.T();
        Algebra::Matrix<double> f1({{1, -1, -2}}), f2({{-1, 2, 3}}), f3({{-1, -1, -1}});
        Algebra::Matrix<double> r1({{-1, -1, -1}}), r2({{4, 3, 1}}), r3({{-3, -2, -1}});
        std::cout << "-\n";
        Algebra::Matrix<double> B({{-2, 5, -6},
                                   {-3, 6, -6},
                                   {-2, 4, -4}});
        std::cout << B * r1.T() << "\n" << B * r2.T() << "\n" << B * r3.T() << "\n";
        Algebra::Matrix<double> A1({{-2, 3,  -5},
                                    {-6, 3,  -9},
                                    {6,  -1, 7}}), B1({{3, 1, 2},
                                                       {3, 0, 3},
                                                       {2, 0, 2}});
        Algebra::Matrix<double> C = A1 * B1;
        std::cout << "C: \n" << C;
        std::cout << "A' : \n" << A1.row_echelon_form() << "\n";
        std::cout << "B' : \n" << B1.row_echelon_form() << "\n";
        //std::cout << A1 * (Algebra::Matrix<double> ({{-1, 1, 1}})).T();
    }

    void Task4() {
        std::cout << "-----Task 4-----\n";
        Rational t = 1;
        std::vector<std::vector<Rational> > arr({{-5,    0,     0, 0, 0},
                                               {1,     -5,    0, 0, 0},
                                               {t - (Rational)1, 0,     4, 0, 1},
                                               {0,     0,     0, 4, 0},
                                               {0,     t + (Rational)1, 0, t, 4}});
        Algebra::Matrix<Rational> A(arr);
        Algebra::Matrix<Rational> A4(A - 4);
        Algebra::Matrix<Rational> A5(A - (-5));
        std::cout << (A5 * A5).rank() << "\n";
        std::cout << (A5 * A5 * A5).rank() << "\n";
        std::cout << (A4).rank() << "\n";
    }

    void Task5() {
        std::cout << "-----Task 5-----\n";
        Algebra::Matrix<double> B({{0, 0, 2,  3},
                                   {0, 0, 3,  5},
                                   {2, 3, 2,  -1},
                                   {3, 5, -1, -7}});
        Algebra::Vector<double> v1({2, -2, 1, -1}), v2({2, 0, 0, 1}), v3({-8, 5, 0, 1}), v4({3, -1, -1, 1});
        std::vector<Algebra::Vector<double> > V({v1, v2, v3, v4});
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                std::cout << V[i].T() * B * V[j] << " ";
            }
            std::cout << "\n";
        }
        Algebra::Matrix<double> Base = v1 | v2 | v3 | v4;
        std::cout << Base.rank();
    }
}

int main() {
    Tasks::Task1();
    Tasks::Task2();
    Tasks::Task3();
    Tasks::Task4();
    Tasks::Task5();
}
