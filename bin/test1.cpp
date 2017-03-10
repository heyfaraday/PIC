#include <iostream>
#include <vector>

#include "boris_push.hpp"
#include "constants.hpp"
#include "field.hpp"
#include "parameters.hpp"

int main() {

    double t = t_start_input;

    std::vector <double> u = U_start;
    std::vector <double> x = X_start;

    std::vector <double> B = B_start;
    std::vector <double> E = E_start;

    FILE *fid = fopen("out.txt", "w");

    while(t < t_end_input)
    {
        fprintf(fid, "%f %f %f %f %f %f %f\n", t, u.at(1), u.at(2), u.at(3), x.at(1), x.at(2), x.at(3));
        
        borisPusher_explicit(u, E, B, chargeOverMass_input, dt_input);
        updateX(u, x, dt_input);
        updateE(E, x, u, t, dt_input);
        updateB(B, x, u, t, dt_input);
        
        t += dt_input;
    }

    return 0;
}