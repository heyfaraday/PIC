#include <vector>
#include "math.h"

#include "parameters.hpp"
#include "constants.hpp"

#include "field.hpp"

void updateE(std::vector<double>& E, std::vector<double> x, std::vector<double> v, double t, double dt) {
    E.at(1) = E.at(1);
    E.at(2) = E.at(2);
    E.at(3) = E.at(3);
}

void updateB(std::vector<double>& B, std::vector<double> x, std::vector<double> v, double t, double dt) {
    B.at(1) = B.at(1);
    B.at(2) = B.at(2);
    B.at(3) = B.at(3);
}