#include "doubledouble.h"
#include <iostream>

int main(int argc, char const *argv[])
{
    DoubleDouble z(100);
    auto r = z.exp();
    std::cout << r << "\n";
    return 0;
}
