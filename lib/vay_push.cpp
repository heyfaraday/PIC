#include <vector>
#include <iostream>
#include "math.h"

#include "utils.hpp"

#include "vay_push.hpp"

void vayPusher(std::vector<double>& v, std::vector<double> E,
               std::vector<double> B, double chargeOverMass, double dt) {
    std::vector<double> u(4);
    std::vector<double> v_cross_B(4);
    std::vector<double> u_cross_tau(4);
    std::vector<double> tau(4);

    double u_star;
    double sigma;

    u.at(1) = v.at(1) * v.at(0);
    u.at(2) = v.at(2) * v.at(0);
    u.at(3) = v.at(3) * v.at(0);

    vectorCross(v, B, v_cross_B);

    u.at(1) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(1) + v_cross_B.at(1));
    u.at(2) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(2) + v_cross_B.at(2));
    u.at(3) += (dt * chargeOverMass / 2.0) * (2.0 * E.at(3) + v_cross_B.at(3));

    tau.at(1) = (dt * chargeOverMass / 2.0) * B.at(1);
    tau.at(2) = (dt * chargeOverMass / 2.0) * B.at(2);
    tau.at(3) = (dt * chargeOverMass / 2.0) * B.at(3);

    u_star = (u.at(1) * tau.at(1) + u.at(2) * tau.at(2) + u.at(3) * tau.at(3));
    sigma = 1 + (u.at(1) * u.at(1) + u.at(2) * u.at(2) + u.at(3) * u.at(3)) -
            (tau.at(1) * tau.at(1) + tau.at(2) * tau.at(2) + tau.at(3) * tau.at(3));

    v.at(0) = sqrt((sigma + sqrt(sigma * sigma + 4 * (tau.at(1) * tau.at(1) + tau.at(2) * tau.at(2) +
            tau.at(3) * tau.at(3) + u_star * u_star))) / 2.0);

    tau.at(1) = tau.at(1) / v.at(0);
    tau.at(2) = tau.at(2) / v.at(0);
    tau.at(3) = tau.at(3) / v.at(0);

    vectorCross(u, tau, u_cross_tau);

    u.at(1) = (1.0 / (1.0 + tau.at(1) * tau.at(1) + tau.at(2) * tau.at(2) + tau.at(3) * tau.at(3)))
              * (u.at(1) + tau.at(1) * (u.at(1) * tau.at(1) + u.at(2) * tau.at(2) + u.at(3) * tau.at(3))
              + u_cross_tau.at(1));

    u.at(2) = (1.0 / (1.0 + tau.at(1) * tau.at(1) + tau.at(2) * tau.at(2) + tau.at(3) * tau.at(3)))
              * (u.at(2) + tau.at(2) * (u.at(1) * tau.at(1) + u.at(2) * tau.at(2) + u.at(3) * tau.at(3))
              + u_cross_tau.at(2));

    u.at(3) = (1.0 / (1.0 + tau.at(1) * tau.at(1) + tau.at(2) * tau.at(2) + tau.at(3) * tau.at(3)))
              * (u.at(3) + tau.at(3) * (u.at(1) * tau.at(1) + u.at(2) * tau.at(2) + u.at(3) * tau.at(3))
              + u_cross_tau.at(3));

    v.at(1) = u.at(1) / v.at(0);
    v.at(2) = u.at(2) / v.at(0);
    v.at(3) = u.at(3) / v.at(0);
}