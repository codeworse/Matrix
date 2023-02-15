#include "polynom.h"
#include<algorithm>

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
    private:
        value *M;
        size_t n, m;

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

        void fail(const char *message) const {
            std::cerr << "FAIL: ";
            std::cerr << message;
            exit(0);
        }

        Matrix<value> E(size_t n_) {
            Matrix<value> ans(n_, n_);
            for (size_t i = 0; i < n_; ++i) {
                value *ans_string = ans.M + i * n_;
                for (size_t j = 0; j < n_; ++j) {
                    *(ans_string + j) = 0;
                }
                *(ans_string + i) = 1;
            }
            return ans;
        }

    public:
        Matrix(const size_t n_ = 1, const size_t m_ = 1, const value &x = value()) : n(n_), m(m_) {
            M = new value[n * m];
            for (size_t i = 0; i < n; ++i) {
                value *now_string = (M + i * m);
                for (size_t j = 0; j < m; ++j) {
                    value &now = *(now_string + j);
                    now = x;
                }
            }
        }

        Matrix(const std::vector<std::vector<value> > &arr) : n(arr.size()), m(arr[0].size()) {
            M = new value[n * m];
            for (size_t i = 0; i < n; ++i) {
                value *now_string = M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(now_string + j) = arr[i][j];
                }
            }
        }

        operator std::vector<std::vector<value> >() const {
            std::vector<std::vector<value> > ans(n, std::vector<value>(m));
            for (size_t i = 0; i < n; ++i) {
                value *now_string = M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    ans[i][j] = *(now_string + j);
                }
            }
            return ans;
        }

        Matrix(const Matrix<value> &other) : n(other.n), m(other.m) {
            M = new value[n * m];
            for (size_t i = 0; i < n; ++i) {
                value *now_string = M + i * m;
                const value *other_string = other.M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(now_string + j) = *(other_string + j);
                }
            }
        }

        Matrix<value> &operator=(const Matrix<value> &other) {
            if (this == &other) {
                return *this;
            }
            delete[] M;
            n = other.n;
            m = other.m;
            M = new value[n * m];
            for (size_t i = 0; i < n; ++i) {
                value *now_string = M + i * m;
                value *other_string = other.M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(M + i * m + j) = *(other_string + j);
                }
            }
            return *this;
        }

        Matrix<value> operator+(const Matrix<value> &other) const {
            if (n != other.n || m != other.m) {
                fail("You cannot add matrix with different sizes");
            }
            Matrix<value> ans(other);
            for (size_t i = 0; i < n; ++i) {
                value *this_string = M + i * m;
                value *ans_string = ans.M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(ans_string + j) += *(this_string + j);
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
                value *other_string = M + i * m;
                value *ans_string = ans.M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(ans_string + j) -= *(other_string + j);
                }
            }
            return ans;
        }

        Matrix<value> operator-() const {
            Matrix<value> ans(*this);
            for (size_t i = 0; i < n; ++i) {
                value *this_string = M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(this_string + j) = -*(this_string + j);
                }
            }
            return ans;
        }

        Matrix<value> operator*(const Matrix<value> &other) const {
            if (m != other.n) {
                fail("You cannot multiply matrix with different sizes");
            }
            Matrix<value> ans(n, other.m);
            for (size_t i = 0; i < n; ++i) {
                value *ans_string = ans.M + i * other.m;
                for (size_t r = 0; r < m; ++r) {
                    value *other_string = other.M + r * other.m;
                    value this_value = *(M + i * m + r);
                    for (size_t j = 0; j < other.m; ++j) {
                        ans_string[j] += this_value * other_string[j];
                    }
                }
            }
            return ans;
        }

        bool operator==(const Matrix<value> &other) const {
            if (n != other.n || m != other.m) {
                return false;
            }
            for (size_t i = 0; i < n; ++i) {
                value *this_string = M + i * m;
                value *other_string = M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    if (!(*(this_string + j) == *(other_string + j))) {
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
                    if (*(ans.M + index * m + i + d) != 0) {
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
                    value *i_string = (ans.M + i * m);
                    value *non_zero_string = (ans.M + non_zero_index * m);
                    for (size_t j = 0; j < m; ++j) {
                        std::swap(*(i_string + j), *(non_zero_string + j));
                    }
                }
                for (size_t j = i + 1; j < n; ++j) {
                    value k = (-(*(ans.M + j * m + i + d))) / *(ans.M + i * m + i + d);
                    for (size_t r = 0; r < m; ++r) {
                        value now_value = *(ans.M + i * m + r);
                        *(ans.M + j * m + r) += k * now_value;
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
            value ans = *(ref.M);
            for (size_t i = 1; i < n; ++i) {
                ans *= *(ref.M + i * m + i);
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
                value *a_string = a.M + i * a.m;
                for (size_t j = 0; j < a.m; ++j) {
                    out << *(a_string + j) << "\t";
                }
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
                        p.get_k(0) = *(M + t[0]);
                        p.get_k(1) = -1;
                    } else {
                        p.get_k(0) = *(M + t[0]);
                    }
                    now = p;
                }
                for (size_t i = 1; i < t.size(); ++i) {
                    Polynom<value> p(2);
                    if (i == t[i]) {
                        p.get_k(0) = *(M + i * m + t[i]);
                        p.get_k(1) = -1;
                    } else {
                        p.get_k(0) = *(M + i * m + t[i]);
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
                value *now_string = Mat.M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    if (*(now_string + j) != 0) {
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
                *(A.M + i * m + i) -= x;
            }
            return A;
        }

        value &get_value(size_t i, size_t j) {
            if (i >= n || j >= m) {
                fail("index out of range in get_value function");
            }
            return *(M + i * m + j);
        }

        Matrix<value> T() const {
            Matrix<value> ans(m, n);
            for (size_t i = 0; i < n; ++i) {
                value *this_string = M + i * m;
                for (size_t j = 0; j < m; ++j) {
                    *(ans.M + j * m + i) = *(this_string + j);
                }
            }
            return ans;
        }

        ~Matrix() {
            delete[] M;
        }
    };


}