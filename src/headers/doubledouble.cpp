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

std::pair<double, double> two_product(const double &x, const double &y)
{
    double u = x * 134217729.0;
    double v = y * 134217729.0;
    double s = u - (u - x);
    double t = v - (v - y);
    double f = x - s;
    double g = y - t;
    double r = x * y;
    double e = ((s * t - r) + s * g + f * t) + f * g;
    return {r, e};
}

DoubleDouble operator+(const DoubleDouble &first, const DoubleDouble &second)
{
    auto [res, err] = two_sum(first.res, second.err);
    err += first.err + second.err;
    auto [new_res, new_err] = two_sum_quick(res, err);
    return DoubleDouble(new_res, new_err);
}

DoubleDouble operator-(const DoubleDouble &first, const DoubleDouble &second)
{
    auto [res, err] = two_differece(first.res, second.res);
    err += first.err + second.res;
    auto [new_res, new_err] = two_sum_quick(res, err);
    return DoubleDouble(new_res, new_err);
}

DoubleDouble operator*(const DoubleDouble &first, const DoubleDouble &second)
{
    auto [res, err] = two_product(first.res, second.res);
    err += first.err + second.err;
    auto [new_res, new_err] = two_sum_quick(res, err);
    return DoubleDouble(new_res, new_err);
}

DoubleDouble operator/(const DoubleDouble &first, const DoubleDouble &second)
{
    double res = first.res / second.res;
    auto [s, f] = two_product(res, second.res);
    double err = (first.res - s - f + first.err - res * second.err) / second.res;
    auto [new_res, new_err] = two_sum_quick(res, err);
    return DoubleDouble(new_res, new_err);
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

// DoubleDouble DoubleDouble::operator+(double other)
// {
//     return *this + DoubleDouble(other);
// }

// DoubleDouble DoubleDouble::operator-(double other)
// {
//     return *this - DoubleDouble(other);
// }

// DoubleDouble DoubleDouble::operator*(double other)
// {
//     return *this * DoubleDouble(other);
// }

// DoubleDouble DoubleDouble::operator/(double other)
// {
//     return *this / DoubleDouble(other);
// }

// DoubleDouble operator+(double left, const DoubleDouble &right)
// {
//     return DoubleDouble(left) + right;
// }

// DoubleDouble operator-(double left, const DoubleDouble &right)
// {
//     return DoubleDouble(left) - right;
// }

// DoubleDouble operator*(double left, const DoubleDouble &right)
// {
//     return DoubleDouble(left) * right;
// }

// DoubleDouble operator/(double left, const DoubleDouble &right)
// {
//     return DoubleDouble(left) / right;
// }

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
DoubleDouble DoubleDouble::round()
{
    return DoubleDouble(std::round(res), std::round(err));
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
    DoubleDouble n = round();
    DoubleDouble x = *this - n;
    DoubleDouble u = (((((((((((x +
                                DoubleDouble(156)) *
                                   x +
                               DoubleDouble(12012)) *
                                  x +
                              DoubleDouble(600600)) *
                                 x +
                             DoubleDouble(21621600)) *
                                x +
                            DoubleDouble(588107520)) *
                               x +
                           DoubleDouble(12350257920)) *
                              x +
                          DoubleDouble(201132771840)) *
                             x +
                         DoubleDouble(2514159648000)) *
                            x +
                        DoubleDouble(23465490048000)) *
                           x +
                       DoubleDouble(154872234316800)) *
                          x +
                      DoubleDouble(647647525324800)) *
                         x +
                     DoubleDouble(1295295050649600);
    DoubleDouble v = (((((((((((x -
                                DoubleDouble(156)) *
                                   x +
                               DoubleDouble(12012)) *
                                  x -
                              DoubleDouble(600600)) *
                                 x +
                             DoubleDouble(21621600)) *
                                x -
                            DoubleDouble(588107520)) *
                               x +
                           DoubleDouble(12350257920)) *
                              x -
                          DoubleDouble(201132771840)) *
                             x +
                         DoubleDouble(2514159648000)) *
                            x -
                        DoubleDouble(23465490048000)) *
                           x +
                       DoubleDouble(154872234316800)) *
                          x -
                      DoubleDouble(647647525324800)) *
                         x +
                     DoubleDouble(1295295050649600);
    return DOUBLEDOUBLE_E.power(n) * (u / v);
}

DoubleDouble::operator double()
{
    return res;
}