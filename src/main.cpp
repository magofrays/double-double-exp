#include "doubledouble.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

int main(int argc, char const *argv[])
{
    std::vector<double> xs;
    std::vector<double> ys;
    json j;
    DoubleDouble s(372.8);
    DoubleDouble e = exp(s);
    std::cout << (std::exp(372.8) - double(e)) << " r: " << e << " e: " << e.err << "\n";
    for (double i = 1; i != 5000; i++)
    {
        double x = i / 10;
        DoubleDouble y(x);
        xs.push_back(x);
        double res = std::abs(std::exp(x) - double(exp(y)));
        if (res > std::pow(10.0, 143))
        {
            std::cout << x << "\n";
        }
        ys.push_back(res);
    }
    j["xs"] = xs;
    j["ys"] = ys;
    std::ofstream file("../test_methods.json");
    if (file.is_open())
    {
        file << j.dump(4);
        file.close();
        std::cout << "JSON saved to ../test_methods.json" << std::endl;
    }
}
