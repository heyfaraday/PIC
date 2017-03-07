#pragma once

void borisPusher(std::vector<double>& u, std::vector<double>& E, std::vector<double>& B,
                 double chargeOverMass, double dt);

void borisPusher_explicit(std::vector<double>& u, std::vector<double>& E, std::vector<double>& B,
                          double charge_over_mass, double dt);

void updateX(std::vector<double>& u, std::vector<double>& x, double dt);

void vectorCross(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);