#include<vector>

template<typename value>
class Polynom {
private:
    std::vector<value> arr;
public:
    Polynom(size_t sz = 0) {
        arr.resize(sz);
    }

    Polynom(const std::initializer_list<value> a) {
        arr.resize(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            arr[i] = a[i];
        }
    }

    Polynom(const Polynom &p) : arr(p.arr) {
    }

    Polynom<value> &operator+=(const Polynom<value> &other) {
        if (arr.size() < other.arr.size()) {
            arr.resize(other.arr.size());
        }
        for (size_t i = 0; i < other.arr.size(); ++i) {
            arr[i] += other.arr[i];
        }
        return *this;
    }

    Polynom<value> &operator-=(const Polynom<value> &other) {
        if (arr.size() < other.arr.size()) {
            arr.resize(other.arr.size());
        }
        for (size_t i = 0; i < arr.size(); ++i) {
            arr[i] -= other.arr[i];
        }
        return *this;
    }

    Polynom<value> operator+(const Polynom<value> &other) {
        Polynom<value> ans(*this);
        ans += other;
        return ans;
    }

    value &get_k(const size_t i) {
        return arr[i];
    }

    Polynom<value> operator-(const Polynom<value> &other) {
        Polynom<value> ans(*this);
        ans -= other;
        return ans;
    }

    friend Polynom<value> operator*(const Polynom<value> &a, const Polynom<value> &b) {
        Polynom<value> ans(a.arr.size() + b.arr.size());
        for (size_t i = 0; i < a.arr.size(); ++i) {
            for (size_t j = 0; j < b.arr.size(); ++j) {
                ans.arr[i + j] += a.arr[i] * b.arr[j];
            }
        }
        return ans;
    }

    Polynom<value> operator*=(const Polynom<value> &other) {
        *this = *this * other;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &out, const Polynom<value> &p) {
        bool first = true;
        for (size_t i = 0; i < p.arr.size(); ++i) {
            if (p.arr[i] == 0)
                continue;
            if (i == 0) {
                out << p.arr[i] << " ";
                first = false;
                continue;
            }
            if (p.arr[i] == 1) {
                if (!first)
                    out << "+ ";
                out << "x^" << i << " ";
                first = false;
                continue;
            }
            if (p.arr[i] == -1) {
                if (first) {
                    out << "-";
                } else {
                    out << "- ";
                }
                out << "x^" << i << " ";
                first = false;
                continue;
            }
            if (p.arr[i] < 0) {
                if (!first)
                    out << "- ";
                else
                    out << "-";
            } else if (!first) {
                out << "+ ";
            }
            out << std::abs(p.arr[i]) << " * x^" << i << " ";
            first = false;
        }
        return out;
    }

};