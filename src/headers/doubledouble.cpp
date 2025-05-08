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

std::pair<double, double> split(const double &a)
{
    double t = 134217729 * a;
    double ahi = t - (t - a);
    double alo = a - ahi;
    return {ahi, alo};
}

std::pair<double, double> two_product(const double &x, const double &y)
{
    double res = x * y;
    auto [xhi, xlo] = split(x);
    auto [yhi, ylo] = split(y);
    double err = ((xhi * yhi - res) + xhi * ylo + xlo * yhi) + xlo * ylo;
    return {res, err};
}

DoubleDouble operator+(const DoubleDouble &first, const DoubleDouble &second)
{
    auto Double = two_sum(first.res, second.res);
    Double.second += first.err + second.err;
    Double = two_sum_quick(Double.first, Double.second);
    return DoubleDouble(Double.first, Double.second);
}

DoubleDouble operator-(const DoubleDouble &first, const DoubleDouble &second)
{
    auto Double = two_differece(first.res, second.res);
    Double.second += (first.err - second.err);
    Double = two_sum_quick(Double.first, Double.second);
    return DoubleDouble(Double.first, Double.second);
}

DoubleDouble operator*(const DoubleDouble &first, const DoubleDouble &second)
{
    auto Double = two_product(first.res, second.res);
    Double.second += first.res * second.err + second.res * first.err;
    Double = two_sum_quick(Double.first, Double.second);
    return DoubleDouble(Double.first, Double.second);
}

DoubleDouble operator/(const DoubleDouble &first, const DoubleDouble &second)
{
    double res = first.res / second.res;
    auto temp = two_product(res, second.res);
    double err = (first.res - temp.first - temp.second + first.err - res * second.err) / second.res;
    auto Double = two_sum_quick(res, err);
    return DoubleDouble(Double.first, Double.second);
}

DoubleDouble &DoubleDouble::operator+=(const DoubleDouble &other)
{
    *this = (static_cast<DoubleDouble>(*this) + other);
    return *this;
}

DoubleDouble &DoubleDouble::operator-=(const DoubleDouble &other)
{
    *this = (static_cast<DoubleDouble>(*this) - other);
    return *this;
}

DoubleDouble &DoubleDouble::operator*=(const DoubleDouble &other)
{
    *this = (static_cast<DoubleDouble>(*this) * other);
    return *this;
}

DoubleDouble &DoubleDouble::operator/=(const DoubleDouble &other)
{
    *this = (static_cast<DoubleDouble>(*this) / other);
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

DoubleDouble DoubleDouble::exp()
{
    DoubleDouble m = ((*this) / DOUBLEDOUBLE_LN2).trunc();
    DoubleDouble x = (*this - DOUBLEDOUBLE_LN2 * m);
    DoubleDouble first(2);
    first = first.power(m);
    DoubleDouble second = DoubleDouble(1) + x + x.power(2) / DoubleDouble(2) + x.power(3) / DoubleDouble(6) +
                          x.power(4) / DoubleDouble(24) + x.power(5) / DoubleDouble(120) + x.power(6) / DoubleDouble(720) +
                          x.power(7) / DoubleDouble(5040) + x.power(8) / DoubleDouble(40320) + x.power(9) / DoubleDouble(362880) + x.power(10) / DoubleDouble(3628800);
    return first * second;
}

DoubleDouble::operator double()
{
    return res;
}