#include "System.h"


System::System(unsigned long n):
Bead(), rand_seed(n), g(n){
  N_springs = N_beads-1;
  f = vecd(N_steps_s, 0);
  vrel_par = vecd(DIM, 0);
  vrel_perp = vecd(DIM, 0);
  v_rel = vecd(DIM, 0);
  u = vecd(DIM, 0);
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

      if (j != 0){
        p[i][j].m = 0.0005;
      } else{
        p[i][j].m = 0.01;
      }
    }
  }
};


void System::initialize(){
  getValues();

  std::normal_distribution<double> distribution(1.0, 0.1);
  std::normal_distribution<double> distribution2(0.0, 0.1);

  for (size_t i = 0; i < N_lines; i++) {
    p[i][0].q[0] = 0.0;
    p[i][0].q[1] = 1.0;
    p[i][0].q[2] = z0;

    for (size_t j = 1; j < N_beads; j++) {
        double a = distribution(g);
        double b = distribution2(g);
        double c = distribution2(g);
        p[i][j].q[0] = s0*j + s0*a;                       //initialize positions
        p[i][j].q[1] = 1.0 + s0*b;
        p[i][j].q[2] = z0 + s0*c;

      for (size_t k = 0; k < DIM; k++) {
        p[i][j].vel[k]=0;
      }

      p[i][j].n=j;

    }

    p[i][0].l=i;
  }

  getWeight();
};


void System::getDistances(){
  double d = 0;
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads-1; j++) {
      d = 0;
      for (size_t k = 0; k < DIM; k++) {
        d += pow((p[i][j+1].q[k]-p[i][j].q[k]), 2);
      }
      p[i][j].r = sqrt(d);
    }
  }
};


void System::getAngles(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads-1; j++) {
      p[i][j].theta = acos((p[i][j+1].q[2] - p[i][j].q[2]) / p[i][j].r);
      p[i][j].phi = atan2((p[i][j+1].q[1] - p[i][j].q[1]), (p[i][j+1].q[0] - p[i][j].q[0]));
    }
  }
};


void System::getFel(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads-1; j++) {
      p[i][j].F_el[0] = k*(p[i][j+1].q[0] - p[i][j].q[0] - s0*cos(p[i][j].phi)*sin(p[i][j].theta));
      p[i][j].F_el[1] = k*(p[i][j+1].q[1] - p[i][j].q[1] - s0*sin(p[i][j].phi)*sin(p[i][j].theta));
      p[i][j].F_el[2] = k*(p[i][j+1].q[2] - p[i][j].q[2] - s0*cos(p[i][j].theta));
    }

    for (size_t k = 0; k < DIM; k++) {
      p[i][N_beads-1].F_el[k] = 0;
    }

    for (size_t j = 1; j < N_beads; j++) {
      p[i][j].F_el[0] -= k*(p[i][j].q[0] - p[i][j-1].q[0] - s0*cos(p[i][j-1].phi)*sin(p[i][j-1].theta));
      p[i][j].F_el[1] -= k*(p[i][j].q[1] - p[i][j-1].q[1] - s0*sin(p[i][j-1].phi)*sin(p[i][j-1].theta));
      p[i][j].F_el[2] -= k*(p[i][j].q[2] - p[i][j-1].q[2] - s0*cos(p[i][j-1].theta));
    }
  }
};


void System::getScalarProd(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 1; j < N_beads-1; j++) {
      p[i][j].sp = 0;
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].sp += (p[i][j-1].q[k] - p[i][j].q[k]) * (p[i][j+1].q[k]-p[i][j].q[k]);
      }
    }
  }
};


void System::getFKP(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 2; j < N_beads; j++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].F_KP[k] = -J*((p[i][j-2].q[k] - p[i][j-1].q[k])/(p[i][j-2].r*(p[i][j-1].r)) - p[i][j-1].sp*(p[i][j].q[k] - p[i][j-1].q[k])/(p[i][j-2].r*pow(p[i][j-1].r, 3)));
      }
    }

    for (size_t k = 0; k < DIM; k++) {
      p[i][0].F_KP[k] = 0;
      p[i][1].F_KP[k] = 0;
    }

    for (size_t j = 1; j < N_beads-1; j++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].F_KP[k] -= J*((2*p[i][j].q[k] - p[i][j-1].q[k] - p[i][j+1].q[k])/(p[i][j-1].r*p[i][j].r) + p[i][j].sp*((p[i][j-1].q[k] - p[i][j].q[k])/(pow(p[i][j-1].r, 3)*p[i][j].r) + (p[i][j+1].q[k] - p[i][j].q[k])/(p[i][j-1].r*pow(p[i][j].r, 3))));
      }
    }

    for (size_t j = 0; j < N_beads-2; j++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].F_KP[k] -= J*((p[i][j+2].q[k] - p[i][j+1].q[k])/(p[i][j].r*p[i][j+1].r) - p[i][j+1].sp*(p[i][j].q[k] - p[i][j+1].q[k])/(pow(p[i][j].r, 3)*(p[i][j+1].r)));
      }
    }
  }
};


void System::getWindVel(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].wind_vel[0] = w*p[i][j].q[2] - h*psi(p[i][j].q[0])*p[i][j].q[2];
      //p[i][j].wind_vel[0] = 1.0;
      p[i][j].wind_vel[1] = 0.0;
      p[i][j].wind_vel[2] = h*d_psi(p[i][j].q[0])*pow(p[i][j].q[2], 2)/2;
      //p[i][j].wind_vel[2] = 0.0;
    }
  }
};


void System::getWeight(){
  p[0][0].F_weight[2] = -p[0][0].m*grav;
};


void System::getDrag(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      for (size_t k = 0; k < DIM; k++) {
        p[i][j].F_drag[k]=-b*(p[i][j].vel[k]-p[i][j].wind_vel[k]);
      }
    }
  }
};


void System::getNewDrag(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 1; j < N_beads-1; j++) {
      norm = 0;
      vrel_par_sp = 0;
      for (size_t k = 0; k < DIM; k++) {
        u[k] = p[i][j+1].q[k] - p[i][j-1].q[k];
        norm += u[k]*u[k];
      }
      for (size_t k = 0; k < DIM; k++) {
        u[k] = u[k]/sqrt(norm);
        v_rel[k] = p[i][j].vel[k] - p[i][j].wind_vel[k];
        vrel_par_sp += v_rel[k]*u[k];
      }
      for (size_t k = 0; k < DIM; k++) {
        vrel_par[k] = vrel_par_sp*u[k];
        vrel_perp[k] = v_rel[k] - vrel_par[k];
        p[i][j].F_drag[k] = -b*vrel_par[k] - 2*b*vrel_perp[k];
      }
    }
  }
};


void System::thermalize(){

  for (size_t t = 0; t < N_steps_therm_i; t++) {
    if (t%frac_dt == 0){
      t_s = t/frac_dt;
      h=f[t_s];
    }

    computeForces();

    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 1; j < N_beads; j++) {
        p[i][j].oneStepProp();
      }
    }

    std::cout << p[2][2].q[0] << '\n';
  }
};


void System::computeForces(){
  getDistances();
  getScalarProd();
  getAngles();
  getWindVel();
  getNewDrag();
  getFKP();
  getFel();
};


void System::evolve(){
  std::ofstream out_x;
  out_x.open("../data/vel_cont0.txt", std::ios::app);
  bool flag = 0;
  initialize();

  thermalize();

  auto start1 = high_resolution_clock::now();
  for (size_t t = N_steps_therm_i; t < N_steps_i - N_steps_therm_i; t++) {

    computeForces();

    //if ((t%frac_dt == 0) && (t < 2*N_steps_therm_i)){
    if (t < 2*N_steps_therm_i){
      out_x << p[0][0].q[0] << " " << p[0][0].vel[0] << "\n";
    }

    if (t%frac_dt == 0){
      t_s = t/frac_dt;
      h=f[t_s];
    }

    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 0; j < N_beads; j++) {

        if ((j == 0) && (i > 0)){
          for (size_t k = 0; k < DIM; k++) {
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
        p[i][0].q[k] = p[N_lines-1][0].q[k];
      }
    }

    if (flag == 1){
      break;
    }
  }

  auto stop1 = high_resolution_clock::now();
  auto duration1 = duration_cast<microseconds>(stop1 - start1);
  std::cout << "Time taken by cycle evolve: " << duration1.count() << " microseconds" << std::endl;


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


double System::d_psi(double x){
  return exp((-pow(x,2)) / (2*pow(sigma_x,2)));
  //return 0;
};


double System::psi(double x){
  return cost_pi*(erf(x/cost)-1);
  //return 0;
};
