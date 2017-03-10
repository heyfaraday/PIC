#import <vector>

#import "math.h"
#import "parameters.hpp"
#import "constants.hpp"

void updateE(std::vector<double>& E, std::vector<double>& x, std::vector<double>& u, double t, double dt) {

    E.at(1) = E.at(1);
    E.at(2) = E.at(2);
    E.at(3) = E.at(3);

}

void updateB(std::vector<double>& B, std::vector<double>& x, std::vector<double>& u, double t, double dt) {

    B.at(1) = B.at(1);
    B.at(2) = B.at(2);
    // B.at(3) = 20.0 - 0.005 * fabs(x.at(1) + x.at(2)); Drift
    B.at(3) = B.at(3);

}