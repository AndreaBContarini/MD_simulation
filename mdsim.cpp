#include "./sim.hpp"
int main(int argc, char **argv)
{
  mdsimLJ<double> md;
  md.init_rng(0); // inizializzo il generatore di numeri casuali
  md.prepare_initial_conf(); // creo la configurazione iniziale
  md.run(); // lancio la simulazione MC
  //md.save_final_results(); // Stampo pressione e Energia finali
  return 0;
}
