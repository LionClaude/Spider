#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <chrono>

using vecd = std::vector<double>;
using vecdp = std::vector<double*>;
using namespace std::chrono;

const double N_sec = 200;
const double N_sec_therm = 10;
const double dt_i = 0.001;
const double dt_s = 0.1;
const int N_steps_i = int (N_sec/dt_i);
const int N_steps_therm_i = int (N_sec_therm/dt_i);
const int N_steps_s = int (N_sec/dt_s);
const int N_steps_therm_s = int (N_sec_therm/dt_s);
const int frac_dt = int (dt_s/dt_i);
const int N_rep = 1;
const int N_beads = 15;
const int N_lines = 3;
const int DIM = 3;
const double k = 20.0;
const double J = 0.01;
const double s0 = 0.1;
const double w = 0.8;
const double grav = 9.81;
const double b = 0.05;
const double z0 = 2.0;
const double sigma_x = 1.0;
const double cost = sigma_x*sqrt(2);
const double cost_pi = sigma_x*sqrt(M_PI/2);
const double lambda = 0.4;
const double sigma = 1.0;

#endif
