#ifndef DOUBLEDOUBLE_H
#define DOUBLEDOUBLE_H
#include <limits>
#include <cmath>

#define DOUBLEDOUBLE_E DoubleDouble(2.718281828459045, 1.4456468917292502e-16)
#define DOUBLEDOUBLE_PI DoubleDouble(3.141592653589793, 1.2246467991473532e-16)
#define DOUBLEDOUBLE_LN2 DoubleDouble(0.6931471805599453, 2.3190468138462996e-17)
#define DOUBLEDOUBLE_INF DoubleDouble(std::numeric_limits<double>::infinity())
#define SPLIT_THRESH 6.69692879491417e+299
#define SPLITTER 134217729.0
#define EPSILON 4.93038065763132e-32
class DoubleDouble
{

public:
    double res, err;
    explicit DoubleDouble(double res, double err = 0.0) : res(res), err(err) {}

    DoubleDouble &operator+=(const DoubleDouble &other);
    DoubleDouble &operator-=(const DoubleDouble &other);
    DoubleDouble &operator*=(const DoubleDouble &other);
    DoubleDouble &operator/=(const DoubleDouble &other);
    DoubleDouble &operator+=(double other);
    DoubleDouble &operator-=(double other);
    DoubleDouble &operator*=(double other);
    DoubleDouble &operator/=(double other);
    operator bool() const;
    bool operator==(const DoubleDouble &other) const;
    bool operator>(const DoubleDouble &other) const;
    bool operator<(const DoubleDouble &other) const;
    bool operator!=(const DoubleDouble &other) const;
    bool operator<=(const DoubleDouble &other) const;
    bool operator>=(const DoubleDouble &other) const;
    operator double();
    DoubleDouble trunc();
    DoubleDouble power(int n);
};

inline DoubleDouble add(double first, double second);
inline DoubleDouble operator+(const DoubleDouble &first, double second);
inline DoubleDouble operator+(double first, const DoubleDouble &second);
inline DoubleDouble sub(double first, double second);
inline DoubleDouble operator-(const DoubleDouble &first, double second);
inline DoubleDouble operator-(double first, const DoubleDouble &second);
inline DoubleDouble mul(double first, double second);
inline DoubleDouble operator*(const DoubleDouble &first, double second);
inline DoubleDouble operator*(double first, const DoubleDouble &second);
inline DoubleDouble div(double first, double second);
inline DoubleDouble operator/(const DoubleDouble &first, double second);
inline DoubleDouble operator/(double first, const DoubleDouble &second);
inline DoubleDouble operator+(const DoubleDouble &first, const DoubleDouble &second);
inline DoubleDouble operator-(const DoubleDouble &first, const DoubleDouble &second);
inline DoubleDouble operator*(const DoubleDouble &first, const DoubleDouble &second);
inline DoubleDouble operator/(const DoubleDouble &first, const DoubleDouble &second);
DoubleDouble exp(const DoubleDouble &a);

static const double inv_fact[15][2] = {
    {1.66666666666666657e-01, 9.25185853854297066e-18},
    {4.16666666666666644e-02, 2.31296463463574266e-18},
    {8.33333333333333322e-03, 1.15648231731787138e-19},
    {1.38888888888888894e-03, -5.30054395437357706e-20},
    {1.98412698412698413e-04, 1.72095582934207053e-22},
    {2.48015873015873016e-05, 2.15119478667758816e-23},
    {2.75573192239858925e-06, -1.85839327404647208e-22},
    {2.75573192239858883e-07, 2.37677146222502973e-23},
    {2.50521083854417202e-08, -1.44881407093591197e-24},
    {2.08767569878681002e-09, -1.20734505911325997e-25},
    {1.60590438368216133e-10, 1.25852945887520981e-26},
    {1.14707455977297245e-11, 2.06555127528307454e-28},
    {7.64716373181981641e-13, 7.03872877733453001e-30},
    {4.77947733238738525e-14, 4.39920548583408126e-31},
    {2.81145725434552060e-15, 1.65088427308614326e-31}};

#endif