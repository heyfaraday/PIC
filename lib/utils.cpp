#include <vector>
#include "math.h"

#include "utils.hpp"

std::vector<double> vectorCross(std::vector<double> a, std::vector<double> b) {

    std::vector<double> result(4);

    result.at(1) = a.at(2) * b.at(3) - a.at(3) * b.at(2);
    result.at(2) = a.at(3) * b.at(1) - a.at(1) * b.at(3);
    result.at(3) = a.at(1) * b.at(2) - a.at(2) * b.at(1);

    return result;
}

double scalarCross(std::vector<double> a, std::vector<double> b) {
    return a.at(1) * b.at(1) + a.at(2) * b.at(2) + a.at(3) * b.at(3);
}

double scalarCross_explicit(std::vector<double> a, std::vector<double> b) {
    return a.at(0) * b.at(0) - a.at(1) * b.at(1) - a.at(2) * b.at(2) - a.at(3) * b.at(3);
}


void updateX(std::vector<double> v, std::vector<double>& x, double dt) {
    x.at(1) += v.at(1) * dt;
    x.at(2) += v.at(2) * dt;
    x.at(3) += v.at(3) * dt;
}

void updateGamma(std::vector<double> v, double& gamma) {
    gamma = 1.0 / sqrt(1 - v.at(1) * v.at(1) - v.at(2) * v.at(2) - v.at(3) * v.at(3));
}

void updateGamma_v(std::vector<double>& v) {
    v.at(0) = 1.0 / sqrt(1 - v.at(1) * v.at(1) - v.at(2) * v.at(2) - v.at(3) * v.at(3));
}