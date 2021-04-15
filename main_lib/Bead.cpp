#include "Bead.h"

Bead::Bead(){
  q_old = vecd(DIM);
  q = vecd(DIM);
  q_new = vecd(DIM);
  vel = vecd(DIM);
  wind_vel = vecd(DIM);
  qlleft = vecdp(DIM);
  qleft = vecdp(DIM);
  qright = vecdp(DIM);
  qrright = vecdp(DIM);
  F_el = vecd(DIM);
  F_KP = vecd(DIM);
  F_drag = vecd(DIM);
  F_weight = vecd(DIM);
  sigma_x = 1.0;
  cost = sigma_x*sqrt(2);
  cost_pi = sigma_x*sqrt(M_PI/2);
};


void Bead::getDistances(){
  r = 0;
  for (size_t i = 0; i < DIM; i++) {
    r += pow((*qright[i] - q[i]), 2);
  }
  r = sqrt(r);
};


void Bead::getAngles(){
  theta = acos((*qright[2] - q[2]) / r);
  phi = atan2((*qright[1] - q[1]), (*qright[0] - q[0]));
};


void Bead::getFel(){
  F_el = vecd(DIM, 0);

  if (n < N_beads - 1){
    F_el[0] = k*(*qright[0] - q[0] - s0*cos(phi)*sin(theta));
    F_el[1] = k*(*qright[1] - q[1] - s0*sin(phi)*sin(theta));
    F_el[2] = k*(*qright[2] - q[2] - s0*cos(theta));
  }
  if (n > 0){
    F_el[0] -= k*(q[0] - *qleft[0] - s0*cos(*phi_neigh)*sin(*theta_neigh));
    F_el[1] -= k*(q[1] - *qleft[1] - s0*sin(*phi_neigh)*sin(*theta_neigh));
    F_el[2] -= k*(q[2] - *qleft[2] - s0*cos(*theta_neigh));
  }
};


void Bead::getScalarProd(){
  sp = 0;

  if ((n > 0) && (n < N_beads-1)){
    for (size_t i = 0; i < DIM; i++) {
      sp += (*qleft[i] - q[i]) * (*qright[i]-q[i]);
    }
  }
};


void Bead::getFKP(){
  F_KP = vecd(DIM, 0);

  for (size_t i = 0; i < DIM; i++) {
    if (n > 1){
      F_KP[i] = -J*((*qlleft[i] - *qleft[i])/(*r_fleft*(*r_left)) - *sp_left*(q[i] - *qleft[i])/(*r_fleft*pow(*r_left, 3)));
    }
    if ((n > 0) && (n < N_beads-1)){
      F_KP[i] -= J*((2*q[i] - *qleft[i] - *qright[i])/(*r_left*r) + sp*((*qleft[i] - q[i])/(pow(*r_left, 3)*r) + (*qright[i] - q[i])/(*r_left*pow(r, 3))));
    }
    if (n < N_beads-2){
      F_KP[i] -= J*((*qrright[i] - *qright[i])/(r*(*r_right)) - *sp_right*(q[i] - *qright[i])/(pow(r, 3)*(*r_right)));
    }
  }
};


void Bead::getVel(){
  for (size_t i = 0; i < DIM; i++) {
    vel[i] = (q_new[i] - q_old[i])/(2*dt_i);
  }
};


void Bead::getWindVel(){
  wind_vel[0] = w*q[2] - h*psi(q[0])*q[2];
  wind_vel[1] = 0.0;
  wind_vel[2] = h*d_psi(q[0])*pow(q[2], 2)/2;
};


void Bead::getWeight(){
  for (size_t i = 0; i < DIM; i++) {
    if ((i == 2) && (n == 0)){
      F_weight[i] = -m*g;
    } else{
      F_weight[i] = 0;
    }
  }
};


void Bead::getDrag(){
  for (size_t i = 0; i < DIM; i++) {
    F_drag[i]=-b*(vel[i]-wind_vel[i]);
  }
};


void Bead::firstStepProp(){
  for (size_t i = 0; i < DIM; i++) {
    q[i] = q_old[i] + vel[i]*dt_i + (pow(dt_i, 2)/(2*m))*(b*(wind_vel[i] - vel[i]) + F_el[i] + F_KP[i] + F_weight[i]);
  }
};


void Bead::oneStepProp(){
  for (size_t i = 0; i < DIM; i++) {
    q_new[i] = (1/(m + b/2*dt_i))*(2*m*q[i] + (b/2*dt_i - m)*q_old[i] + pow(dt_i, 2)*(F_el[i] + b*wind_vel[i] + F_KP[i] + F_weight[i]));
    q_old[i] = q[i];
    q[i] = q_new[i];
  }
};


double Bead::d_psi(double x){
  //return exp((-pow(x,2)) / (2*pow(sigma_x,2)));
  return 0;
};


double Bead::psi(double x){
  //return cost_pi*(erf(x/cost)-1);
  return 0;
};
