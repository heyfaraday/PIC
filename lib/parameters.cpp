#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "parameters.hpp"

bool is_number(const std::string& s) {

    if (s.empty())
        return false;
    bool trigger = false;
    if (s[0] == '.') {
        return false;
    }
    for (int i = 0; i < s.size(); i++) {
        char c = s[i];
        if (c == '.' && !trigger) {
            trigger = !trigger;
            continue;
        } else if (c == '.') {
            return false;
        }
        if (i == 0 && c == '-') {
            continue;
        }
        if (!std::isdigit(c))
            return false;
    }
    return true;
}

double read_from_file(const std::string param) {

    std::ifstream infile("bin/input.parameters");
    std::string str;

    if (infile.is_open()) {

        while (getline(infile, str)) {

            std::istringstream iss(str);
            std::vector<std::string> words{std::istream_iterator<std::string>{iss},
                                           std::istream_iterator<std::string>{}};
            if (words.size() < 1)
                continue;
            if (words[0] == param) {
                if (is_number(words[1])) {
                    return atof(words[1].c_str());
                } else {
                    std::cout << "Error converting string to number for parameter '"
                              << param << "' in the input file.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
    } else {
        std::cout << "Cannot open the file.\n";
        exit(EXIT_FAILURE);
    }
    std::cout << "Cannot find the given parameter '" << param
              << "' in the input file.\n";
    exit(EXIT_FAILURE);
}

const double Bx = read_from_file("Bx_input");
const double By = read_from_file("By_input");
const double Bz = read_from_file("Bz_input");

const double Ex = read_from_file("Ex_input");
const double Ey = read_from_file("Ey_input");
const double Ez = read_from_file("Ez_input");

const double X1 = read_from_file("X1_input");
const double X2 = read_from_file("X2_input");
const double X3 = read_from_file("X3_input");

const double U1 = read_from_file("U1_input");
const double U2 = read_from_file("U2_input");
const double U3 = read_from_file("U3_input");

const double dt_input = read_from_file("dt_input");
const double t_start_input = read_from_file("t_start_input");
const double t_end_input = read_from_file("t_end_input");
const double chargeOverMass_input = read_from_file("chargeOverMass_input");

const std::vector<double> B_start{0.0, Bx, By, Bz};
const std::vector<double> E_start{0.0, Ex, Ey, Ez};

const std::vector<double> X_start{0.0, X1, X2, X3};
const std::vector<double> U_start{0.0, U1, U2, U3};