#include "doubledouble.h"
#include <utility>
#include <cmath>

std::pair<double, double> two_sum_quick(const double &x, const double &y)
{
    double r = x + y;
    double e = y - (r - x);
    return {r, e};
}

std::pair<double, double> two_sum(const double &x, const double &y)
{
    double r = x + y;
    double t = r - x;
    double e = (x - (r - t)) + (y - t);
    return {r, e};
}

std::pair<double, double> two_differece(const double &x, const double &y)
{
    double r = x - y;
    double t = r - x;
    double e = (x - (r - t)) - (y + t);
    return {r, e};
}

std::pair<double, double> split(double a)
{
    double temp;
    double hi, lo;
    if (a > SPLIT_THRESH || a < -SPLIT_THRESH)
    {
        a *= 3.7252902984619140625e-09; // 2^-28
        temp = SPLITTER * a;
        hi = temp - (temp - a);
        lo = a - hi;
        hi *= 268435456.0; // 2^28
        lo *= 268435456.0; // 2^28
    }
    else
    {
        temp = SPLITTER * a;
        hi = temp - (temp - a);
        lo = a - hi;
    }
    return {hi, lo};
}

std::pair<double, double> two_product(const double &x, const double &y)
{
    double res = x * y;
    auto [xhi, xlo] = split(x);
    auto [yhi, ylo] = split(y);
    double err = ((xhi * yhi - res) + xhi * ylo + xlo * yhi) + xlo * ylo;
    return {res, err};
}

// ADDITION
inline DoubleDouble add(double first, double second)
{
    auto res = two_sum(first, second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator+(const DoubleDouble &first, double second)
{
    auto res = two_sum(first.res, second);
    res.second += first.err;
    res = two_sum_quick(res.first, res.second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator+(double first, const DoubleDouble &second)
{
    return (second + first);
}

inline DoubleDouble operator+(const DoubleDouble &first, const DoubleDouble &second) // IEEE ADD
{
    auto Res = two_sum(first.res, second.res);
    auto Error = two_sum(first.err, second.err);
    Res.second += Error.first;
    Res = two_sum_quick(Res.first, Res.second);
    Res.second += Error.second;
    Res = two_sum_quick(Res.first, Res.second);
    return DoubleDouble(Res.first, Res.second);
}

// SUBTRACT
inline DoubleDouble sub(double first, double second)
{
    auto res = two_differece(first, second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator-(const DoubleDouble &first, double second)
{
    auto res = two_differece(first.res, second);
    res.second += first.err;
    res = two_sum_quick(res.first, res.second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator-(double first, const DoubleDouble &second)
{
    auto res = two_differece(first, second.res);
    res.second += second.err;
    res = two_sum_quick(res.first, res.second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator-(const DoubleDouble &first, const DoubleDouble &second) // IEEE SUB
{
    auto Res = two_differece(first.res, second.res);
    auto Err = two_differece(first.err, second.err);
    Res.second += Err.first;
    Res = two_sum_quick(Res.first, Res.second);
    Res.second += Err.second;
    Res = two_sum_quick(Res.first, Res.second);
    return DoubleDouble(Res.first, Res.second);
}

// MULTIPLY
inline DoubleDouble mul(double first, double second)
{
    auto res = two_product(first, second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator*(const DoubleDouble &first, double second)
{
    auto res = two_product(first.res, second);
    res.second += (first.err * second);
    res = two_sum_quick(res.first, res.second);
    return DoubleDouble(res.first, res.second);
}

inline DoubleDouble operator*(double first, const DoubleDouble &second)
{
    return (second * first);
}

inline DoubleDouble operator*(const DoubleDouble &first, const DoubleDouble &second)
{
    auto Double = two_product(first.res, second.res);
    Double.second += first.res * second.err + second.res * first.err;
    Double = two_sum_quick(Double.first, Double.second);
    return DoubleDouble(Double.first, Double.second);
}

// DIVISION
inline DoubleDouble div(double a, double b)
{

    double q1 = a / b;

    /* Compute  a - q1 * b */
    auto res = two_product(q1, b);
    auto s = two_differece(a, res.first);
    s.second -= res.second;

    /* get next approximation */
    double q2 = (s.first + s.second) / b;

    s = two_sum_quick(q1, q2);

    return DoubleDouble(s.first, s.second);
}

inline DoubleDouble operator/(const DoubleDouble &first, double second)
{
    double q1 = first.res / second;
    auto p = two_product(q1, second);
    auto s = two_differece(first.res, p.first);
    s.second += first.err;
    s.second -= p.second;
    double q2 = (s.first + s.second) / second;
    auto r = two_sum_quick(q1, q2);
    return DoubleDouble(r.first, r.second);
}

inline DoubleDouble operator/(double first, const DoubleDouble &second)
{
    return DoubleDouble(first) / second;
}

inline DoubleDouble operator/(const DoubleDouble &first, const DoubleDouble &second)
{
    double q1 = first.res / second.res;
    DoubleDouble r = first - q1 * second;
    double q2 = r.res / second.res;
    r -= (q2 * second);
    double q3 = r.res / second.res;
    auto res = two_sum_quick(q1, q2);
    r = DoubleDouble(res.first, res.second) + q3;
    return r;
}

DoubleDouble &DoubleDouble::operator+=(const DoubleDouble &other)
{
    *this = *this + other;
    return *this;
}

DoubleDouble &DoubleDouble::operator-=(const DoubleDouble &other)
{
    *this = *this - other;
    return *this;
}

DoubleDouble &DoubleDouble::operator*=(const DoubleDouble &other)
{
    *this = *this * other;
    return *this;
}

DoubleDouble &DoubleDouble::operator/=(const DoubleDouble &other)
{
    *this = *this / other;
    return *this;
}

DoubleDouble &DoubleDouble::operator+=(double other)
{
    *this = *this + other;
    return *this;
}

DoubleDouble &DoubleDouble::operator-=(double other)
{
    *this = *this - other;
    return *this;
}

DoubleDouble &DoubleDouble::operator*=(double other)
{
    *this = *this * other;
    return *this;
}

DoubleDouble &DoubleDouble::operator/=(double other)
{
    *this = *this / other;
    return *this;
}

DoubleDouble::operator bool() const
{
    return res != 0.0 && err != 0.0;
}

bool DoubleDouble::operator==(const DoubleDouble &other) const
{
    return res == other.res && err == other.err;
}

bool DoubleDouble::operator!=(const DoubleDouble &other) const
{
    return res != other.res || err != other.err;
}

bool DoubleDouble::operator>(const DoubleDouble &other) const
{
    return res > other.res || res == other.res && err > other.err;
}

bool DoubleDouble::operator>=(const DoubleDouble &other) const
{
    return res > other.res || res == other.res && err >= other.err;
}

bool DoubleDouble::operator<(const DoubleDouble &other) const
{
    return res < other.res || res == other.res && err < other.err;
}

bool DoubleDouble::operator<=(const DoubleDouble &other) const
{
    return res < other.res || res == other.res && err <= other.err;
}
DoubleDouble DoubleDouble::trunc()
{
    return DoubleDouble(std::trunc(res), std::trunc(err));
}

DoubleDouble DoubleDouble::power(int n)
{
    /* Пример
    Шаг 1: i=1101 → res *= b (b^1)
    Шаг 2: i=110  → b = b^2
    Шаг 3: i=11   → res *= b^2 (теперь res = b^1 * b^2 = b^3)
    Шаг 4: i=1    → res *= b^4 (итого res = b^3 * b^4 = b^7)
    */
    DoubleDouble b = *this;
    int i = std::abs(n);
    DoubleDouble res(1.0);
    while (true)
    {
        if ((i & 1) == 1)
            res *= b;
        if (i <= 1)
            break;
        i >>= 1;
        b *= b;
    }
    if (n < 0)
        return DoubleDouble(1.0) / res;
    return res;
}

inline DoubleDouble mul_pwr2(const DoubleDouble &a, double b) // double-double * double,  where double is a power of 2.
{
    return DoubleDouble(a.res * b, a.err * b);
}

inline std::pair<double, double> two_sqr(double a) // Computes fl(a*a) and err(a*a)
{
    double res = a * a;
    auto [hi, lo] = split(a);
    double err = ((hi * hi - res) + 2.0 * hi * lo) + lo * lo;
    return {res, err};
}

inline DoubleDouble sqr(const DoubleDouble &a) // Squaring
{
    auto p = two_sqr(a.res);
    p.second += 2.0 * a.res * a.err;
    p.second += a.err * a.err;
    auto s = two_sum_quick(p.first, p.second);
    return DoubleDouble(s.first, s.second);
}

inline DoubleDouble ldexp(const DoubleDouble &a, int exp)
{
    return DoubleDouble(std::ldexp(a.res, exp), std::ldexp(a.err, exp));
}

DoubleDouble exp(const DoubleDouble &a)
{
    double k = 512.0;
    double inv_k = 1.0 / k;
    if (a.res <= -709.0)
    {
        return DoubleDouble(0.0);
    }
    if (a.res >= 709.0)
    {
        return DOUBLEDOUBLE_INF;
    }
    if (a.res == 0.0)
    {
        return DoubleDouble(1.0);
    }
    if (a.res == 1.0)
    {
        return DOUBLEDOUBLE_E;
    }
    double m = std::floor(a.res / DOUBLEDOUBLE_LN2 + 0.5);
    DoubleDouble r = mul_pwr2(a - DOUBLEDOUBLE_LN2 * m, inv_k);
    DoubleDouble p = sqr(r);
    DoubleDouble s = r + mul_pwr2(p, 0.5);
    p *= r;
    DoubleDouble t = p * DoubleDouble(inv_fact[0][0], inv_fact[0][1]);
    int i = 0;
    do
    {
        s += t;
        p *= r;
        ++i;
        t = p * DoubleDouble(inv_fact[i][0], inv_fact[i][1]);
    } while (std::abs(double(t)) > inv_k * EPSILON && i < 5);
    s += t;
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s += 1.0;
    return ldexp(s, static_cast<int>(m));
}

DoubleDouble::operator double()
{
    return res;
}