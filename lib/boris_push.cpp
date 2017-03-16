#include <vector>
#include "math.h"

#include "utils.hpp"

#include "boris_push.hpp"

void borisPusher(std::vector<double>& v, std::vector<double>& E,
                 std::vector<double>& B, double chargeOverMass, double dt) {
    std::vector<double> u_s(4);
    std::vector<double> u_cross_b(4);
    double f1, f2;

    v.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    v.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    v.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);

    f1 = (dt * chargeOverMass) / (2.0);
    f2 = (2.0 * f1) /
         (1.0 + f1 * f1 * (B.at(1) * B.at(1) + B.at(2) * B.at(2) + B.at(3) * B.at(3)));

    vectorCross(v, B, u_cross_b);

    u_s.at(1) = v.at(1) + f1 * u_cross_b.at(1);
    u_s.at(2) = v.at(2) + f1 * u_cross_b.at(2);
    u_s.at(3) = v.at(3) + f1 * u_cross_b.at(3);

    vectorCross(u_s, B, u_cross_b);

    v.at(1) += f2 * u_cross_b.at(1);
    v.at(2) += f2 * u_cross_b.at(2);
    v.at(3) += f2 * u_cross_b.at(3);

    v.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    v.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    v.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);
}

void borisPusher_explicit(std::vector<double>& v, std::vector<double>& E,
                          std::vector<double>& B, double chargeOverMass, double dt) {
    std::vector<double> u_s(4);
    std::vector<double> u_cross_b(4);
    double f1, f2;

    v.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    v.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    v.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);

    v.at(0) = sqrt(1.0 + v.at(1) * v.at(1) + v.at(2) * v.at(2) + v.at(3) * v.at(3));

    f1 = (dt * chargeOverMass) / (2.0 * v.at(0));
    f2 = (2.0 * f1) /
         (1.0 + f1 * f1 * (B.at(1) * B.at(1) + B.at(2) * B.at(2) + B.at(3) * B.at(3)));

    vectorCross(v, B, u_cross_b);

    u_s.at(1) = v.at(1) + f1 * u_cross_b.at(1);
    u_s.at(2) = v.at(2) + f1 * u_cross_b.at(2);
    u_s.at(3) = v.at(3) + f1 * u_cross_b.at(3);

    vectorCross(u_s, B, u_cross_b);

    v.at(1) += f2 * u_cross_b.at(1);
    v.at(2) += f2 * u_cross_b.at(2);
    v.at(3) += f2 * u_cross_b.at(3);

    v.at(0) = 0.0;
    v.at(1) += (dt * chargeOverMass / 2.0) * E.at(1);
    v.at(2) += (dt * chargeOverMass / 2.0) * E.at(2);
    v.at(3) += (dt * chargeOverMass / 2.0) * E.at(3);
}