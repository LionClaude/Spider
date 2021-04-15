#include "System.h"


System::System(unsigned long n):
Bead(), rand_seed(n), g(n){
  N_springs = N_beads-1;
  z0 = 2.0;
  lambda = 0.4;
  sigma = 1.0;
  f = vecd(N_steps_s, 0);
};


void System::genNumber(){
  std::poisson_distribution<int> poisson(lambda*N_sec);
  N_events = poisson(g);
};


void System::genPoint(){
  std::uniform_real_distribution<double> uniform(0.0,N_sec);
  mu = uniform(g);
};


void System::genGaussian(){
  std::normal_distribution<double> gauss(0.0,1.0);
  double a = gauss(g);
  for (size_t i = 0; i < N_steps_s; i++) {
    t = i*dt_s;
    f[i] += a * exp(-pow(t-mu,2) / (2*pow(sigma,2)));
  }
};


void System::genSeries(){
  genNumber();
  for (size_t i = 0; i < N_events; i++) {
    genPoint();
    genGaussian();
  }
  printSeries();
};


void System::printSeries(){
  std::ofstream out_w;
  out_w.open("../data/windvel_start0.txt", std::ios::app);

  for (size_t t = 0; t < N_steps_therm_s; t++) {
    out_w << f[t] << "\n";
  }
  out_w << "\n";
  out_w.close();
};


void System::getValues(){

  p = std::vector<std::vector<Bead> > (N_lines, std::vector<Bead> (N_beads));

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].h = f[0];

      if (j != 0){
        p[i][j].m = 0.0005;
      } else{
        p[i][j].m = 0.01;
      }

      if ((j == 0) && (i > 0)){
        p[i][0].g=0.0;
        p[i][0].b=0.0;
      } else{
        p[i][j].b=0.1;
        p[i][j].g=9.81;
      }
    }
  }
};


void System::initialize(){
  getValues();

  std::normal_distribution<double> distribution(1.0, 0.1);
  std::normal_distribution<double> distribution2(0.0, 0.1);

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      if ((j == 0) && (i > 0)){
        for (size_t k = 0; k < DIM; k++) {
          p[i][j].q_old[k] = p[0][0].q_old[k];
        }
      }
      else{
        double a = distribution(g);
        double b = distribution2(g);
        double c = distribution2(g);
        p[i][j].q_old[0] = s0*j + s0*a;                       //initialize positions
        p[i][j].q_old[1] = 1.0 + s0*b;
        p[i][j].q_old[2] = z0 + s0*c;
      }

      for (size_t k = 0; k < DIM; k++) {
        p[i][j].q[k]=p[i][j].q_old[k];
        p[i][j].vel[k]=0;
      }
      p[i][j].n=j;
    }
    p[i][0].l=i;
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].qleft[k]=&p[i][(j-1+N_beads)%N_beads].q[k];       //initialize pointers to neighbors with PBC
        p[i][j].qlleft[k]=&p[i][(j-2+N_beads)%N_beads].q[k];
        p[i][j].qright[k]=&p[i][(j+1)%N_beads].q[k];
        p[i][j].qrright[k]=&p[i][(j+2)%N_beads].q[k];
      }
      p[i][j].getDistances();
      p[i][j].getScalarProd();
    }
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].r_left=&p[i][(j-1+N_beads)%N_beads].r;
      p[i][j].r_right=&p[i][(j+1)%N_beads].r;
      p[i][j].r_fleft=&p[i][(j-2)%N_beads].r;
      p[i][j].sp_left=&p[i][(j-1+N_beads)%N_beads].sp;
      p[i][j].sp_right=&p[i][(j+1)%N_beads].sp;
      p[i][j].getAngles();
    }
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].theta_neigh=&p[i][(j-1+N_beads)%N_beads].theta;
      p[i][j].phi_neigh=&p[i][(j-1+N_beads)%N_beads].phi;
    }
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].getWindVel();
      p[i][j].getWeight();
      p[i][j].getFKP();
      p[i][j].getFel();
    }
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      if ((j == 0) && (i > 0)){
        for (size_t k = 0; k < DIM; k++) {
          p[i][0].q_old[k]=p[i-1][0].q[k];
        }
        p[i][j].firstStepProp();
      }
    }
  }

  for (size_t i = 0; i < N_lines - 1; i++) {
    for (size_t k = 0; k < DIM; k++) {
      p[i][0].q_old[k]=p[N_lines-1][0].q_old[k];
      p[i][0].q[k]=p[N_lines-1][0].q[k];
    }
  }
};


void System::thermalize(){

  for (size_t t = 0; t < N_steps_therm_i; t++) {
    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 0; j < N_beads; j++) {
        if (t%frac_dt == 0){
          t_s = t/frac_dt;
          p[i][j].h=f[t_s];
        }
      }
    }

    computeForces();

    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 1; j < N_beads; j++) {
        p[i][j].oneStepProp();
      }
    }
  }
};


void System::computeForces(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].getDistances();
      p[i][j].getScalarProd();
      p[i][j].getAngles();
      p[i][j].getWindVel();
      p[i][j].getVel();
      p[i][j].getDrag();
      p[i][j].getWeight();
      //p[i][j].getFKP();
      p[i][j].getFel();
    }
  }
};


void System::evolve(){
  std::ofstream out_x;
  out_x.open("../data/vel_cont0.txt", std::ios::app);
  bool flag = 0;
  initialize();
  thermalize();

  for (size_t t = N_steps_therm_i; t < N_steps_i - N_steps_therm_i; t++) {
    computeForces();
    if ((t%frac_dt == 0) && (t < 2*N_steps_therm_i)){
    //if (t < 5*N_steps_therm_i){
      out_x << p[0][0].q[0] << " " << p[0][0].q_old[0] << " " << (p[0][0].q[0]-p[0][0].q_old[0])/dt_i << "\n";
    }
    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 0; j < N_beads; j++) {
        if (t%frac_dt == 0){
          t_s = t/frac_dt;
          p[i][j].h=f[t_s];
        }
        if ((j == 0) && (i > 0)){
          for (size_t k = 0; k < DIM; k++) {
            p[i][0].q_old[k] = p[i - 1][0].q_old[k];
            p[i][0].q[k] = p[i - 1][0].q[k];
          }
        }

        p[i][j].oneStepProp();

        if (p[i][j].q[2] < 0.0){
          flag = 1;
        }
      }
    }

    for (size_t i = 0; i < N_lines - 1; i++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][0].q_old[k] = p[N_lines-1][0].q_old[k];
        p[i][0].q[k] = p[N_lines-1][0].q[k];
      }
    }

    if (flag == 1){
      break;
    }

    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 0; j < N_beads; j++) {
        p[i][j].getVel();
      }
    }
  }

  //printCoord();
  out_x << "\n";
  out_x.close();

};


void System::printCoord(){
  std::ofstream outdata;
  outdata.open("../data/coord50.txt");

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      for (size_t k = 0; k < DIM; k++) {
        outdata << p[i][j].q[k] << " ";
      }
      outdata << "\n";
    }
    outdata << "\n\n";
  }

  outdata.close();
};
