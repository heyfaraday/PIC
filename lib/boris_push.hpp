#pragma once

void borisPusher(std::vector<double>& v, std::vector<double> E, std::vector<double> B,
                 double chargeOverMass, double dt);

void borisPusher_explicit(std::vector<double>& v, std::vector<double> E, std::vector<double> B,
                          double charge_over_mass, double dt);