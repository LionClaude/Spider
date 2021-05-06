#include "Bead.h"

Bead::Bead(){
  q = vecd(DIM, 0);
  vel = vecd(DIM, 0);
  wind_vel = vecd(DIM, 0);
  acc = vecd(DIM, 0);
  F_el = vecd(DIM, 0);
  F_KP = vecd(DIM, 0);
  F_drag = vecd(DIM, 0);
  F_weight = vecd(DIM, 0);
  r = 0;
  theta = 0;
  phi = 0;
  sp = 0;
};


void Bead::oneStepProp(){
  if (n == 0){
    for (size_t i = 0; i < DIM; i++) {
      acc[i] = (F_el[i] + F_KP[i] + 1/N_lines*F_weight[i])/m;
      vel[i] += acc[i]*dt_i;
      q[i] += vel[i]*dt_i;
    }
  } else{
    for (size_t i = 0; i < DIM; i++) {
      acc[i] = (F_el[i] + F_KP[i] + F_drag[i])/m;
      vel[i] += acc[i]*dt_i;
      q[i] += vel[i]*dt_i;
    }
  }
};
