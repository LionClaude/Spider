#ifndef SYSTEM_H
#define SYSTEM_H

#include "Bead.h"


class System: public Bead{
  public:
    int N_springs;
    std::vector<std::vector<Bead> > p;
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
    void getDistances();
    void getAngles();
    void getFel();
    void getScalarProd();
    void getFKP();
    void getWindVel();
    void getDrag();
    void getWeight();
    void getNewDrag();

    double psi(double x);
    double d_psi(double x);

    double T_tot;
    int N_events;
    unsigned long n;
    unsigned long rand_seed;
    double mu;
    double h;
    double t;
    double norm;
    double vrel_par_sp;
    int t_s;
    vecd f;
    vecd vrel_par;
    vecd vrel_perp;
    vecd v_rel;
    vecd u;
    std::default_random_engine g;
    System(unsigned long n);
};

#endif
