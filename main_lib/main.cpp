#include "System.h"


int main(int argc, char const *argv[]) {
  std::ofstream outdata;
  outdata.open("../data/distances0.txt");

  std::cout << "Algorithm started" << '\n';
  for (size_t i = 0; i < N_rep; i++) {
    //unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned long seed = 1234;
    System S(seed);
    S.genSeries();
    S.evolve();

    outdata << S.p[0][0].q[0] << "\n";

    /*if (i % (N_rep/10) == 0){
      std::cout << (float(i)/N_rep)*100 << " % completed" << '\n';
    }*/
  }

  outdata.close();
  return 0;
};
