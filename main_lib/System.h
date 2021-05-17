#ifndef SYSTEM_H
#define SYSTEM_H

#include "Bead.h"


class System: public Bead{

  public:
    System(unsigned long n);
    void evolve();
    void genSeries();
    std::vector<std::vector<Bead> > p;

  private:
    int N_springs;
    void getValues();
    void computeForces();
    void initialize();
    void thermalize();
    void printCoord();
    void genNumber();
    void genPoint();
    void genGaussian();
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
    void addNoise();
    double psi(double x);
    double d_psi(double x);

    double T_tot;
    int N_events;
    unsigned long n;
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
    //std::poisson_distribution<int> poisson;
    //std::uniform_real_distribution<double> uniform;
    //std::normal_distribution<double> gauss;
};

#endif
