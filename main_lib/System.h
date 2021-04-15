#ifndef SYSTEM_H
#define SYSTEM_H

#include "Bead.h"


class System: public Bead{
  public:
    int N_springs;
    std::vector<std::vector<Bead> > p;
    double z0;
    void getValues();
    void computeForces();
    void initialize();
    void thermalize();
    void evolve();
    void printCoord();
    void genNumber();
    void genPoint();
    void genGaussian();
    void genSeries();
    void printSeries();
    void printVel();
    double T_tot;
    int N_events;
    unsigned long n;
    unsigned long rand_seed;
    double lambda;
    double mu;
    double sigma;
    double t;
    int t_s;
    vecd f;
    std::default_random_engine g;
    System(unsigned long n);
};

#endif
