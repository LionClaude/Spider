#ifndef BEAD_H
#define BEAD_H

#include "constants.h"


class Bead{
  public:
    vecd q;
    vecd vel;
    vecd acc;
    vecd wind_vel;
    vecd F_el;
    vecd F_KP;
    vecd F_drag;
    vecd F_weight;
    double r;                  //distances between beads
    double sp;                 //scalar products
    double theta;
    double phi;
    double h;
    double m;
    void oneStepProp();
    int i;
    int n;
    int l;
    Bead();
};

#endif
