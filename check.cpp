#include<iostream>
#include<vector>
#include "matrix.h"
#include "rational.h"

namespace Matrix_tests {
    void fail(const char *message) {
        std::cerr << "FAIL: ";
        std::cerr << message;
        exit(0);
    }

    void check_constructor() {
        std::cout << "Start checking constructor..." << std::endl;
        {
            Algebra::Matrix<long long> Mat_long(2, 2, 2);
            std::vector<std::vector<long long> > Mat_long_vec = Mat_long;
            if (std::vector<std::vector<long long> >({{2, 2},
                                                      {2, 2}}) != Mat_long_vec) {
                fail("wrong convert to std::vector<std::vector<long long> >");
            }
        }
        std::cout << "Finish checking constructor" << std::endl;
    }

    void check_multiplication() {
        std::cout << "Start checking multiplication..." << std::endl;
        {
            Algebra::Matrix<int> A({{1, 2},
                                    {4, 3},
                                    {5, 1}});
            Algebra::Matrix<int> B({{0, 4,  3, 10},
                                    {3, 12, 8, 3}});
            Algebra::Matrix<int> correct_ans({{6, 28, 19, 16},
                                              {9, 52, 36, 49},
                                              {3, 32, 23, 53}});
            Algebra::Matrix<int> C = A * B;
            if (!(C == correct_ans)) {
                fail("wrong answer in multiplication");
            }
        }
        std::cout << "Finish checking multiplication" << std::endl;
    }

    void check_REF() {
        std::cout << "Start checking REF..." << std::endl;
        {
            Algebra::Matrix<double> A({{1, 2, 3},
                                       {4, 5, 6},
                                       {7, 8, 10}});
            auto A_REF = A.row_echelon_form();
            std::vector<std::vector<double> > res = A_REF;
            for (auto &l: res) {
                for (auto x: l) {
                    std::cout << x << " ";
                }
                std::cout << "\n";
            }
            std::cout << A_REF.det() << "\n";
        }
        {
            Algebra::Matrix<double> A({{1, 1, 1},
                                       {1, 1, 1},
                                       {1, 1, 1}});
            std::vector<std::vector<double> > res = A.row_echelon_form();
        }
        {

        }
        std::cout << "Finish checking REF" << std::endl;
    }

    void check_det() {
        std::cout << "Start checking det..." << std::endl;
        std::cout << "Finish checking det..." << std::endl;
    }

    void run_all() {
        check_constructor();
        check_multiplication();
        check_REF();
        check_det();
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
        auto u3 = u1.T();
        std::cout << "ok\n";
        Algebra::Matrix<double> R = A - 3;
        Algebra::Matrix<double> ans = R * u3;
        std::cout << ans;
    }

    void Task4() {
        std::cout << "-----Task 4-----\n";

        Algebra::Matrix<double> A({{3,  -4, -1, 4},
                                   {-7, 2,  4,  -3},
                                   {1,  -4, 1,  4},
                                   {-9, 0,  5,  -2}});
        std::cout << "rk(A) = " << A.rank() << "\n";
        std::cout << "rk(A^2) = " << (A * A).rank() << "\n";
        std::cout << "rk(A^3) = " << (A * A * A).rank() << "\n";

        std::cout << "A^2:\n" << A * A;

        std::cout << "rk(A - 2) = " << (A - 2).rank() << "\n";
        std::cout << "rk((A - 2)^2) = " << ((A - 2) * (A - 2)).rank() << "\n";
        std::cout << "rk((A - 2)^3) = " << ((A - 2) * (A - 2) * (A - 2)).rank() << "\n";
        std::cout << "(A - 2)^2:\n";
        std::cout << ((A - 2) * (A - 2)) << "\n";

        Algebra::Matrix<double> B({{1, 1},
                                   {0, -2},
                                   {2, 0},
                                   {0, -2}});
        Algebra::Matrix<double> C({{1, 0, -1, 0},
                                   {0, 4, 0,  -3}});
        std::cout << B * C << "\n";
    }
}

int main() {
    Tasks::Task1();
    Tasks::Task2();
}
