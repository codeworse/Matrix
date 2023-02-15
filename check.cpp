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
            Algebra::Matrix<int> M_int;
            Algebra::Matrix<long long> Mat_long(2, 2, 2);
            std::vector<std::vector<long long> > Mat_long_vec = Mat_long;
            if (std::vector<std::vector<long long>>{{2, 2},
                                                    {2, 2}} != Mat_long_vec) {
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
    void Task2() {
        std::cout << "-----Task 2-----\n";
        Algebra::Matrix<double> A({{-4, 12,  3,  -9},
                                   {-4, 6,   -3, 3},
                                   {12, -18, 9,  -9},
                                   {4,  -6,  3,  -3}});
        std::cout << "A.T -> \n" << (A.T()).row_echelon_form() << "\n";
        std::cout << "A -> \n" << A.row_echelon_form() << "\n";
        Algebra::Matrix<double> B({{-9, -4, 4, 0},
                                   {15, 8,  0, 4},
                                   {-1, 0,  0, 0}});
        std::cout << "rk(B) = " << B.rank() << "\n";
    }

    void Task3() {
        std::cout << "-----Task 3-----\n";

        Algebra::Matrix<double> A({{4, -4, 9,  5},
                                   {2, -2, 11, 6},
                                   {0, 0,  4,  2},
                                   {0, 0,  -4, -2}});
        std::cout << (A * A) << "\n\n";
        std::cout << (A * A).row_echelon_form() << "\n\n";
        std::cout << "rk(A) = " << A.rank() << "\n";
        std::cout << "rk(A^2) = " << (A * A).rank() << "\n";
        std::cout << "rk(A^3) = " << (A * A * A).rank() << "\n";
        std::cout << "rk(A^4) = " << (A * A * A * A).rank() << "\n";
        Algebra::Matrix<double> A_T = A.T();
        std::cout << (A_T * A_T) << "\n\n";
        std::cout << (A_T * A_T).row_echelon_form() << "\n\n";
        std::cout << "rk(A_T) = " << A_T.rank() << "\n";
        std::cout << "rk(A_T^2) = " << (A_T * A_T).rank() << "\n";
        std::cout << "rk(A_T^3) = " << (A_T * A_T * A_T).rank() << "\n";
        std::cout << "rk(A_T^4) = " << (A_T * A_T * A_T * A_T).rank() << "\n";
        Algebra::Matrix<double> C({{4,  -4, 0, 1},
                                   {0,  0,  2, 1},
                                   {1,  -2, 3, 0},
                                   {-1, 2,  0, 3}});
        std::cout << "rk(C) = " << C.rank() << "\n";
        Algebra::Matrix<double> Cv({{4,  -4, 0, 1, -2},
                                    {0,  0,  2, 1, -9},
                                    {1,  -2, 3, 0, -5},
                                    {-1, 2,  0, 3, -7}});
        std::cout << Cv.row_echelon_form() << "\n";
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
    Tasks::Task2();
    Tasks::Task3();
    Tasks::Task4();
}
