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

const int N_sec = 200;
const int N_sec_therm = 10;
const double dt_i = 0.01;
const double dt_s = 0.1;
const int N_steps_i = int (N_sec/dt_i);
const int N_steps_therm_i = int (N_sec_therm/dt_i);
const int N_steps_s = int (N_sec/dt_s);
const int N_steps_therm_s = int (N_sec_therm/dt_s);
const int frac_dt = int (dt_s/dt_i);
const int N_rep = 50;
const int N_beads = 15;
const int N_lines = 3;
const int DIM = 3;
const double k = 0.5;
const double J = 0.0;
const double s0 = 0.1;
const double w = 0.8;

#endif
