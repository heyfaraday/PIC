#import <vector>

#import "math.h"
#import "parameters.hpp"
#import "constants.hpp"

void updateE(std::vector<double>& E, double t, double dt) {

    E.at(1) = E.at(1);
    E.at(2) = sin(0.5*t);
    E.at(3) = E.at(3);

}

void updateB(std::vector<double>& B, double t, double dt) {

    B.at(1) = B.at(1);
    B.at(2) = B.at(2);
    B.at(3) = cos(0.5*t);

}