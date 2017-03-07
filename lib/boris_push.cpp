#import "math.h"
#include <iostream>
#include <vector>

#include "boris_push.hpp"

void vectorCross(std::vector<double>& a, std::vector<double>& b,
                 std::vector<double>& result) {
    result.at(1) = a.at(2) * b.at(3) - a.at(3) * b.at(2);
    result.at(2) = a.at(3) * b.at(1) - a.at(1) * b.at(3);
    result.at(3) = a.at(1) * b.at(2) - a.at(2) * b.at(1);
}

void borisPusher(std::vector<double>& u, std::vector<double>& E,
                 std::vector<double>& B, double chargeOverMass, double dt) {

    std::vector<double> u_s(4);
    std::vector<double> u_cross_b(4);
    double f1, f2;

    u.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    u.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    u.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);

    f1 = (dt * chargeOverMass) / (2.0);
    f2 = (2.0 * f1) /
         (1.0 + f1 * f1 * (B.at(1) * B.at(1) + B.at(2) * B.at(2) + B.at(3) * B.at(3)));

    vectorCross(u, B, u_cross_b);

    u_s.at(1) = u.at(1) + f1 * u_cross_b.at(1);
    u_s.at(2) = u.at(2) + f1 * u_cross_b.at(2);
    u_s.at(3) = u.at(3) + f1 * u_cross_b.at(3);

    vectorCross(u_s, B, u_cross_b);

    u.at(1) += f2 * u_cross_b.at(1);
    u.at(2) += f2 * u_cross_b.at(2);
    u.at(3) += f2 * u_cross_b.at(3);

    u.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    u.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    u.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);
}

void borisPusher_explicit(std::vector<double>& u, std::vector<double>& E,
                          std::vector<double>& B, double chargeOverMass, double dt) {
    std::vector<double> u_s(4);
    std::vector<double> u_cross_b(4);
    double f1, f2;

    u.at(1) += (dt * chargeOverMass / 2.0) * E.at(0);
    u.at(2) += (dt * chargeOverMass / 2.0) * E.at(1);
    u.at(3) += (dt * chargeOverMass / 2.0) * E.at(2);

    u.at(0) = sqrt(1.0 + u.at(1) * u.at(1) + u.at(2) * u.at(2) + u.at(3) * u.at(3));

    f1 = (dt * chargeOverMass) / (2.0 * u.at(0));
    f2 = (2.0 * f1) /
         (1.0 + f1 * f1 * (B.at(1) * B.at(1) + B.at(2) * B.at(2) + B.at(3) * B.at(3)));

    vectorCross(u, B, u_cross_b);

    u_s.at(1) = u.at(1) + f1 * u_cross_b.at(1);
    u_s.at(2) = u.at(2) + f1 * u_cross_b.at(2);
    u_s.at(3) = u.at(3) + f1 * u_cross_b.at(3);

    vectorCross(u_s, B, u_cross_b);

    u.at(1) += f2 * u_cross_b.at(1);
    u.at(2) += f2 * u_cross_b.at(2);
    u.at(3) += f2 * u_cross_b.at(3);

    u.at(0) = 0.0;
    u.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    u.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    u.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);
}

void updateX(std::vector<double>& u, std::vector<double>& x, double dt) {

    x.at(1) += u.at(1) * dt;
    x.at(2) += u.at(2) * dt;
    x.at(3) += u.at(3) * dt;
}