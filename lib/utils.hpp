#pragma once

std::vector<double> vectorCross(std::vector<double> a, std::vector<double> b);

double scalarCross(std::vector<double> a, std::vector<double> b);

double scalarCross_explicit(std::vector<double> a, std::vector<double> b);

void updateX(std::vector<double> v, std::vector<double>& x, double dt);

void updateGamma(std::vector<double> v, double& gamma);

void updateGamma_v(std::vector<double>& v);