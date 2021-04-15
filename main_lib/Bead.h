#ifndef BEAD_H
#define BEAD_H

#include "constants.h"


class Bead{
  public:
    vecd q_old;
    vecd q;
    vecd q_new;
    vecd vel;
    vecd wind_vel;
    vecd F_el;
    vecd F_KP;
    vecd F_drag;
    vecd F_weight;
    double r;                  //distances between beads
    double *r_left;
    double *r_right;
    double *r_fleft;
    double sp;                 //scalar products
    double *sp_left;
    double *sp_right;
    double theta;
    double phi;
    double *theta_neigh;
    double *phi_neigh;
    double h;
    double b;
    double m;
    double g;
    double sigma_x;
    double cost;
    double cost_pi;
    void getDistances();
    void getAngles();
    void getFel();
    void getScalarProd();
    void getFKP();
    void getVel();
    void getWindVel();
    void getDrag();
    void getWeight();
    void oneStepProp();
    void firstStepProp();
    double psi(double x);
    double d_psi(double x);
    int i;
    int n;
    int l;
    vecdp qlleft;
    vecdp qleft;
    vecdp qright;
    vecdp qrright;
    Bead();
};

#endif
