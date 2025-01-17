#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"

template<typename ntype>
class simparsLJ
{
public:
  int nx, ny, nz; /* nx*ny*nz particelle */
  ntype T, P; // temperature and pressure
  int Np; // numero di particelle
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
  ntype rho, rc; // density & cutoff radius rc
  int simtype; // simulation type (see below) to choose
  int seed; // -1 means random
  pvector<ntype, 3> L; // box
  ntype sigma, epsilon; // Lennard-Jones parameters
  ntype dt, deltra, vmax; // parameter of MC moves
  simparsLJ()
    {
      simtype = 0; // 0 NTV, 1 NPT
      nx = 8;  // number of particles along each direction
      ny = 8;
      nz = 8;
      //Np = nx*ny*nz; //non necessario, scritto già in sim.hpp
      sigma=1.0;
      epsilon=1.0;
      rho = 0.4;
      rc = 2.5 * sigma; //raggio di cutoff settato a 2.5 volte sigma
      seed=0;
      adjstps = -200;
      maxadjstps = 2000;
      eqstps = 5000;


      /////////////////////////////
      totsteps = 10000;
      dt = 0.004;
      /////////////////////////////

      save_mgl_snapshot = 1000;
      savemeasure=20;
      outstps=200;
      T = 1.4;
      P = 0.1035; //se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella 
      deltra = 0.2;
      vmax=10.0;
    }
};
#endif
