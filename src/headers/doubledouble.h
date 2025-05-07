#ifndef DOUBLEDOUBLE_H
#define DOUBLEDOUBLE_H

#define DOUBLEDOUBLE_E DoubleDouble(2.718281828459045, 1.4456468917292502e-16)
#define DOUBLEDOUBLE_PI DoubleDouble(3.141592653589793, 1.2246467991473532e-16)

class DoubleDouble
{
    double res, err;

public:
    explicit DoubleDouble(double res, double err = 0.0) : res(res), err(err) {}
    friend DoubleDouble operator+(const DoubleDouble &first, const DoubleDouble &second);
    friend DoubleDouble operator-(const DoubleDouble &first, const DoubleDouble &second);
    friend DoubleDouble operator*(const DoubleDouble &first, const DoubleDouble &second);
    friend DoubleDouble operator/(const DoubleDouble &first, const DoubleDouble &second);
    DoubleDouble &operator+=(const DoubleDouble &other);
    DoubleDouble &operator-=(const DoubleDouble &other);
    DoubleDouble &operator*=(const DoubleDouble &other);
    DoubleDouble &operator/=(const DoubleDouble &other);
    // DoubleDouble operator+(double other);
    // DoubleDouble operator-(double other);
    // DoubleDouble operator*(double other);
    // DoubleDouble operator/(double other);
    // friend DoubleDouble operator+(double left, const DoubleDouble &right);
    // friend DoubleDouble operator-(double left, const DoubleDouble &right);
    // friend DoubleDouble operator*(double left, const DoubleDouble &right);
    // friend DoubleDouble operator/(double left, const DoubleDouble &right);
    operator bool() const;
    bool operator==(const DoubleDouble &other) const;
    bool operator>(const DoubleDouble &other) const;
    bool operator<(const DoubleDouble &other) const;
    bool operator!=(const DoubleDouble &other) const;
    bool operator<=(const DoubleDouble &other) const;
    bool operator>=(const DoubleDouble &other) const;
    operator double();
    DoubleDouble round();
    DoubleDouble power(int n);
    DoubleDouble exp();
};

#endif