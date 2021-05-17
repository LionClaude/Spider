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
  //std::default_random_engine u;
  //std::normal_distribution<double> noise;
  if (n == 0){
    //std::cout << F_weight[2] << " " << F_el[2] << " " << acc[2] << " " << vel[2] << " " << q[2] << '\n';
    for (size_t i = 0; i < DIM; i++) {
      acc[i] = (F_el[i] + F_KP[i] + F_weight[i]/N_lines)/m;
      vel[i] += acc[i]*dt_i;
      q[i] += vel[i]*dt_i;
    }
  } else{
    for (size_t i = 0; i < DIM; i++) {
      acc[i] = (F_el[i] + F_KP[i] + F_drag[i])/m;
      //acc[i] += noise(u);
      //std::cout << noise(u) << '\n';
      vel[i] += acc[i]*dt_i;
      q[i] += vel[i]*dt_i;
    }
  }
};
