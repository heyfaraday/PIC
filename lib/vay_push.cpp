#include <vector>
#include <iostream>
#include "math.h"

#include "utils.hpp"

#include "vay_push.hpp"

void vayPusher(std::vector<double>& v, std::vector<double> E,
               std::vector<double> B, double chargeOverMass, double dt) {
    std::vector<double> u(4);
    std::vector<double> u_cross_tau(4);
    std::vector<double> tau(4);

    double u_star;
    double sigma;

    u.at(1) = v.at(1) * v.at(0);
    u.at(2) = v.at(2) * v.at(0);
    u.at(3) = v.at(3) * v.at(0);

    u.at(1) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(1) + vectorCross(v, B).at(1));
    u.at(2) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(2) + vectorCross(v, B).at(2));
    u.at(3) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(3) + vectorCross(v, B).at(3));

    tau.at(1) = (dt * chargeOverMass / 2.0) * B.at(1);
    tau.at(2) = (dt * chargeOverMass / 2.0) * B.at(2);
    tau.at(3) = (dt * chargeOverMass / 2.0) * B.at(3);

    u_star = scalarCross(u, tau);
    sigma = 1 + scalarCross(u, u) - scalarCross(tau, tau);

    v.at(0) = sqrt((sigma + sqrt(sigma * sigma + 4 * (scalarCross(tau, tau) + u_star * u_star))) / 2.0);

    tau.at(1) = tau.at(1) / v.at(0);
    tau.at(2) = tau.at(2) / v.at(0);
    tau.at(3) = tau.at(3) / v.at(0);

    u.at(1) = (1.0 / (1.0 + scalarCross(tau, tau))) * (u.at(1) + tau.at(1) * (scalarCross(u, tau))
                                                       + vectorCross(u, tau).at(1));

    u.at(2) = (1.0 / (1.0 + scalarCross(tau, tau))) * (u.at(2) + tau.at(2) * (scalarCross(u, tau))
                                                       + vectorCross(u, tau).at(2));

    u.at(3) = (1.0 / (1.0 + scalarCross(tau, tau))) * (u.at(3) + tau.at(3) * (scalarCross(u, tau))
                                                       + vectorCross(u, tau).at(3));

    v.at(1) = u.at(1) / v.at(0);
    v.at(2) = u.at(2) / v.at(0);
    v.at(3) = u.at(3) / v.at(0);

}