#include "polynom.h"
#include<algorithm>
#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H
namespace Algebra {
    template<typename value>
    class Matrix;

    template<typename value>
    Matrix<value> eye(size_t n) {
        Matrix<value> A(n, n);
        for (size_t i = 0; i < n; ++i) {
            A.get_value(i, i) = 1;
        }
        return A;
    }

    template<typename value>
    class Matrix {
    protected:
        std::vector<std::vector<value> > M;
        size_t n = 1, m = 1;

        bool sign(const std::vector<size_t> &t) const {
            int c = 0;
            for (size_t i = 0; i < t.size(); ++i) {
                for (size_t j = i + 1; j < t.size(); ++j) {
                    if (t[i] > t[j]) {
                        c++;
                    }
                }
            }
            return !(c % 2);
        }

        static void fail(const char *message) {
            std::cerr << "FAIL: ";
            std::cerr << message;
            exit(0);
        }

        Matrix<value> E(size_t n_) {
            Matrix<value> ans(n_, n_);
            for (size_t i = 0; i < n_; ++i) {
                for (size_t j = 0; j < n_; ++j) {
                   M[i][j] = 0;
                }
                M[i][i] = 1;
            }
            return ans;
        }

    public:
        Matrix() {
            n = 1;
            m = 1;
            M.resize(n, std::vector<value> (m, value()));
        }
        Matrix(const size_t n_, const size_t m_, const value x = value()) : n(n_), m(m_) {
            M.resize(n, std::vector<value> (m));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    M[i][j] = x;
                }
            }
        }

        Matrix(const std::vector<std::vector<value> > &arr) : n(arr.size()), m(arr[0].size()) {
            M.resize(n, std::vector<value> (m));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    M[i][j] = arr[i][j];
                }
            }
        }
        operator std::vector<std::vector<value> >() const {
            return M;
        }

        Matrix(const Matrix<value> &other) : n(other.n), m(other.m) {
            M.resize(n, std::vector<value> (m));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    M[i][j] = other.M[i][j];
                }
            }
        }

        Matrix<value> &operator=(const Matrix<value> &other) {
            if (this == &other) {
                return *this;
            }
            n = other.n;
            m = other.m;
            M.resize(n, std::vector<value> (m));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    M[i][j] = other.M[i][j];
                }
            }
            return *this;
        }

        Matrix<value> operator+(const Matrix<value> &other) const {
            if (n != other.n || m != other.m) {
                fail("You cannot add matrix with different sizes");
            }
            Matrix<value> ans(*this);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    ans.M[i][j] += other.M[i][j];
                }
            }
            return ans;
        }

        Matrix<value> operator-(const Matrix<value> &other) const {
            if (n != other.n || m != other.m) {
                fail("You cannot subtract matrix with different sizes");
            }
            Matrix<value> ans(*this);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    ans.M[i][j] -= other.M[i][j];
                }
            }
            return ans;
        }

        Matrix<value> operator-() const {
            Matrix<value> ans(*this);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    M[i][j] = -M[i][j];
                }
            }
            return ans;
        }

        Matrix<value> operator*(const Matrix<value> &other) const {
            if (m != other.n) {
                fail("You cannot multiply matrix with different sizes");
            }
            Matrix<value> ans(n, other.m, 0);
            for (size_t i = 0; i < n; ++i) {
                for (size_t r = 0; r < m; ++r) {
                    value this_value = M[i][r];
                    for (size_t j = 0; j < other.m; ++j) {
                        ans.M[i][j] += this_value * other.M[r][j];
                    }
                }
            }
            return ans;
        }
        Matrix<value> operator*(const value &x) {
            Matrix<value> ans(*this);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    ans.M[i][j] *= x;
                }
            }
            return ans;
        }
        bool operator==(const Matrix<value> &other) const {
            if (n != other.n || m != other.m) {
                return false;
            }
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    if (!(M[i][j] == other.M[i][j])) {
                        return false;
                    }
                }
            }
            return true;
        }

        Matrix<value> row_echelon_form() const {
            Matrix<value> ans(*this);
            size_t d = 0;
            for (int i = 0; i < n && i < m; ++i) {
                int non_zero_index = i - 1;
                if (i + d >= m) {
                    break;
                }
                for (int index = i; index < n; ++index) {
                    if (ans.M[index][i + d] != 0) {
                        non_zero_index = index;
                        break;
                    }
                }
                if (non_zero_index < i) {
                    d++;
                    if (d == m) {
                        break;
                    }
                    i--;
                    continue;
                }
                if (non_zero_index != i) {
                    for (size_t j = 0; j < m; ++j) {
                        std::swap(ans.M[i][j], ans.M[non_zero_index][j]);
                    }
                }
                for (size_t j = i + 1; j < n; ++j) {
                    value k = (-ans.M[j][i + d]) / ans.M[i][i + d];
                    for (size_t r = 0; r < m; ++r) {
                        value now_value = ans.M[i][r];
                        ans.M[j][r] += k * now_value;
                    }
                }
            }
            return ans;
        }

        value det() const {
            if (n != m) {
                fail("det can be count only for square matrix");
            }
            Matrix<value> ref = this->row_echelon_form();
            value ans = ref.M[0][0];
            for (size_t i = 1; i < n; ++i) {
                ans *= ref.M[i][i];
            }
            return ans;
        }

        Matrix<value> power(size_t p) {
            if (n != m) {
                fail("Sizes of matrix in power func are not equal");
            }
            if (p == 1) {
                return *this;
            }
            if (p == 0) {
                return E(n);
            }
            Matrix<value> ans = power(p / 2);
            if (p % 2) {
                return ans * ans * (*this);
            }
            return ans * ans;
        }

        friend std::ostream &operator<<(std::ostream &out, const Matrix<value> &a) {
            for (size_t i = 0; i < a.n; ++i) {
                for (size_t j = 0; j < a.m; ++j) {
                    out << a.M[i][j] << "\t";
                }
                if (i < a.n - 1)
                    out << "\n";
            }
            return out;
        }

        Polynom<value> Characteristic_poly() const {
            if (n != m) {
                fail("Not square matrix in Characteristic_poly");
            }
            Polynom<value> poly;
            std::vector<size_t> t(n);
            for (size_t i = 0; i < n; ++i) {
                t[i] = i;
            }
            bool run = true;
            while (run) {
                Polynom<value> now;
                {
                    Polynom<value> p(2);
                    if (0 == t[0]) {
                        p.get_k(0) = M[0][0];
                        p.get_k(1) = -1;
                    } else {
                        p.get_k(0) = M[0][t[0]];
                    }
                    now = p;
                }
                for (size_t i = 1; i < t.size(); ++i) {
                    Polynom<value> p(2);
                    if (i == t[i]) {
                        p.get_k(0) = M[i][t[i]];
                        p.get_k(1) = -1;
                    } else {
                        p.get_k(0) = M[i][t[i]];
                    }
                    now = now * p;
                }
                if (sign(t))
                    poly += now;
                else
                    poly -= now;
                run = std::next_permutation(t.begin(), t.end());
            }
            return poly;
        }

        int rank() const {
            Matrix<value> Mat = this->row_echelon_form();
            int rk = n;
            for (size_t i = 0; i < n; ++i) {
                bool zero_line = true;
                for (size_t j = 0; j < m; ++j) {
                    if (Mat.M[i][j] != 0) {
                        zero_line = false;
                    }
                }
                if (zero_line)
                    rk--;
            }
            return rk;
        }

        Matrix<value> operator-(const value x) const {
            if (n != m) {
                fail("You cannot subtract value from non-square matrix");
            }
            Matrix<value> A(*this);
            for (size_t i = 0; i < n; ++i) {
                A.M[i][i] -= x;
            }
            return A;
        }

        value &get_value(size_t i, size_t j) {
            if (i >= n || j >= m) {
                fail("index out of range in get_value function");
            }
            return M[i][j];
        }

        Matrix<value> T() const {
            Matrix<value> ans(m, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    ans.M[j][i] = M[i][j];
                }
            }
            return ans;
        }
        friend Matrix<value> operator|(const Matrix<value> &a, const Matrix<value> &b) {
            if (a.n != b.n) {
                fail("in operator| diff sizes");
            }
            Matrix<value> ans(a.n, a.m + b.m);
            for (size_t i = 0; i < a.n; ++i) {
                for (size_t j = 0; j < a.m; ++j) {
                    ans.M[i][j] = a.M[i][j];
                }
                for (size_t j = 0; j < b.m; ++j) {
                    ans.M[i][j + a.m] = b.M[i][j];
                }
            }
            return ans;
        }
        size_t get_column() const {
            return n;
        }
        size_t get_row() const {
            return m;
        }
        ~Matrix() {
        }
    };
}
#endif MATRIX_MATRIX_H