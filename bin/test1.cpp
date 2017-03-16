#include <iostream>
#include <vector>

#include "constants.hpp"
#include "parameters.hpp"

#include "boris_push.hpp"
#include "vay_push.hpp"
#include "utils.hpp"
#include "field.hpp"

int main() {

    double t = t_start_input;

    std::vector<double> v = V_start;
    std::vector<double> x = X_start;

    std::vector<double> B = B_start;
    std::vector<double> E = E_start;

    updateGamma_v(v);

    FILE *fid = fopen("out.txt", "w");

    while (t < t_end_input) {
        vayPusher(v, E, B, chargeOverMass_input, dt_input);
        updateX(v, x, dt_input);

        fprintf(fid, "%f %f %f %f %f %f %f\n", t, v.at(1), v.at(2), v.at(3),
                x.at(1), x.at(2), x.at(3));

        t += dt_input;
    }
    return 0;
}