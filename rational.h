#include<iostream>
#ifndef MATRIX_RATIONAL_H
#define MATRIX_RATIONAL_H
class Rational {
private:
    long long x, y; /// x / y
    long long gcd(long long a, long long b) {
        if (!b)
            return a;
        return gcd(b, a % b);
    }

    void update() {
        if (y < 0) {
            x *= -1;
            y *= -1;
        }
        long long c = gcd(std::abs(x), y);
        x /= c;
        y /= c;
    }

public:
    Rational(const long long a = 0, const long long b = 1) : x(a), y(b) {
    }

    Rational(const Rational &other) : x(other.x), y(other.y) {
    }

    Rational &operator=(const Rational &other) {
        if (this == &other) {
            return *this;
        }
        x = other.x;
        y = other.y;
        return *this;
    }

    Rational &operator+=(const Rational &other) {
        x = x * other.y + other.x * y;
        y *= other.y;
        update();
        return *this;
    }

    Rational &operator-=(const Rational &other) {
        x = x * other.y - other.x * y;
        y *= other.y;
        update();
        return *this;
    }

    Rational operator+(const Rational &other) {
        Rational ans(*this);
        ans += other;
        return ans;
    }

    Rational operator-(const Rational &other) {
        Rational ans(*this);
        ans -= other;
        return ans;
    }

    Rational &operator*=(const Rational &other) {
        x *= other.x;
        y *= other.y;
        update();
        return *this;
    }

    Rational operator*(const Rational &other) {
        Rational ans(*this);
        ans *= other;
        return ans;
    }

    Rational &operator/=(const Rational &other) {
        x *= other.y;
        y *= other.x;
        update();
        return *this;
    }

    Rational operator/(const Rational &other) {
        Rational ans(*this);
        ans /= other;
        return ans;
    }

    Rational reverse() const {
        Rational ans(y, x);
        ans.update();
        return ans;
    }

    Rational operator-() const {
        Rational ans(-x, y);
        ans.update();
        return ans;
    }
    friend Rational operator-(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans = ans - b;
        return ans;
    }
    friend Rational operator+(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans = ans + b;
        return ans;
    }
    friend std::istream &operator>>(std::istream &in, Rational &r) {
        in >> r.x >> r.y;
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const Rational &r) {
        out << r.x << "/" << r.y;
        return out;
    }
    operator double() const {
        return static_cast<double>(x) / y;
    }
    long long getx() const {
        return x;
    }
    long long gety() const {
        return y;
    }
};
namespace std {
    Rational abs(const Rational &a) {
        Rational ans(std::abs(a.getx()), a.gety());
        return ans;
    }
}
#endif //MATRIX_RATIONAL_H