#ifndef _SIMCLASS_
#define _SIMCLASS_
/*
 *  base class (common stuff for MC and MD) 
 *   [X]  calcenergyi()
 *   [X]  calctotenergy()
 *   [X]  pbc() [ to apply pbc ]
 *   [X]  prepare_initial_conf() 
 *   [X]  init_rng()
 *   [X]  run() method (main loop) 
 *   [X]  save_mgl_snapshot() save configurations in mgl format
 *
 ** MC class
 *
 *  NTV ensemble
 *
 *   [X]  trial_move()
 *   [X]  acc_move()
 *   [X]  move_NTV()
 
 *  NPT ensemble
 *
 *   [X]  trial_move_box()
 *   [X]  acc_move_box()
 *   [X]  move_box()
 *   [X]  store_all_pars()
 *   [X]  restore_all_pars()
 *  
 *   [X]  calc_acceptance and adjust()
 *   [X]  init_measures() truncate measure files
 *   [X]  save_measures() append measures
 *
 ** MD class
 *
 *   [X]  gauss() (gauss stochastic variable using Box-Muller)
 *   [*]  calc_forces()
 *   [*]  run (velocity verlet as trotter factorization of propagator)
 *   [*]  prepare_initial_conf(void) (initialize velocities)
 *   [X]  init_measures() truncate measure files
 *   [X]  save_measures() append measures
 *
 *  optional exercises:
 *
 *   [ ]  save_snapshot() save configurations <--- mostrare Lab 28/11/23
 *   [ ]  save_restart()
 */
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <iomanip> // for setprecision()
template <typename ntype>
class simLJ
{
  using simp = simparsLJ<ntype>;
  using pLJ = particleLJ<ntype>;
protected:
  simp pars;  // parametri per composizione
  std::vector<pLJ> parts; // vettore di tutte le particelle

  ntype calcenergyi(int i)
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < pars.Np; j++) // pars.Np è il numero totale di particelle
        {
          if (i==j)
            continue;
          enei += parts[i].vij(parts[j], pars.L); 
          // pars.L è un vettore con i lati del box 
        }
      return enei;
    }
  ntype totenergy(void)
    {
      int i;
      ntype ene=0.0;
      for (i=0; i < pars.Np; i++)
        {
          ene += calcenergyi(i);
        }
      return ene*0.5; //avoid double counting
    }
  void pbc(int i)
    {
      auto Dr = parts[i].r; //contiene la posizione della particella i
      Dr = pars.L.mulcw(rint(Dr.divcw(pars.L))); // L*rint(Dr/L) -> è la PBC
      parts[i].r = parts[i].r - Dr; //aggiorno posizione particella sottraendo la correzione Dr
    }

  void save_mgl_snapshot(long int t)
    {
       std::fstream f;
       std::string s;

      s = "cnf-" + std::to_string(t) + ".mgl";

       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << parts[i].r[0] << " " << parts[i].r[1] << " " <<
             parts[i].r[2] << " @ " << pars.sigma*0.5 << "\n";
         }
       f.close();
    }
 
public:
  void set_sim_type(int type)
    {
      pars.simtype = type;
    }
 
  void init_rng(int n)
    {
      if (n < 0)
        rng.rseed();
      else
        rng.seed(n);
    }
  void prepare_initial_conf(void)
    {
      int ix, iy, iz;
      int cc=0;
      pars.Np = pars.nx*pars.ny*pars.nz;
      parts.resize(pars.Np);
      ntype vcell = pow(pars.sigma,3.0);
      ntype rhomax = 1.0/vcell;
      ntype sf;
      sf = cbrt(rhomax/pars.rho);
      pars.L << pars.nx, pars.ny, pars.nz;
      ntype clen = sf*pars.sigma;
      pars.L *= clen;
      for (ix = 0; ix < pars.nx; ix++)
        for (iy = 0; iy < pars.ny; iy++)
          for (iz = 0; iz < pars.nz; iz++)
            {
              parts[cc].r << ix*clen, iy*clen, iz*clen;
              parts[cc].r -= pars.L*0.5;
              parts[cc].set_sigma(pars.sigma);
              parts[cc].set_epsilon(pars.epsilon);
              parts[cc].set_rcut(pars.rc);
              cc++;
            }
          }


  void run(void) 
    {
      // intentionally void
    }; 
};

template <typename ntype>
class mcsimLJ: public simLJ<ntype>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  long int tot_tra, tot_vol, tra_rej, vol_rej;
  using bc=simLJ<ntype>;
  using bc::pars;
  using bc::calcenergyi, bc::totenergy;
  using bc::parts, bc::pbc;
  using bc::save_mgl_snapshot;

  void trial_move(int i)
    {
      // muovere la particella usando come max 
      // displacement pars.deltra
      // applico le pbc (tramite il metodo pbc)
      pvector<ntype,3> delr;
      delr << pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5);
      parts[i].tra_move(delr);
      pbc(i);
    }
  void acc_move(int i, ntype eno)
    {
      // assume kB=1 (Boltzmann constants) so that beta=1/T 
      // calcolo la nuova energia 
      // accetto la mossa con probabilità min{1taU)}
      // se deltaU < 0 accetto la mossa
      // altrimenti genero xi un numero a caso tra 0 e 1 
      // e se xi < exp(-beta*deltaU) accetto altrimenti rifiuto
      // accetta la mossa con criterio Metropolis MC
      ntype enn = calcenergyi(i);
      ntype delu= enn-eno;
      ntype xi = rng.ranf();
      if (delu > 0.0 && xi >= exp(-delu/pars.T))
        {
          // reject move
          tra_rej++;
          parts[i].restore();
        }
    }
  
  void move_NTV(int i)
    {
      // 1) calcolare l'energia della particella i-esima
      // 2) memorizzo la posizione della particella i
      // 3) trial move 
      // 4) acceptance per la particella i
      ntype eno;
      eno = calcenergyi(i);
      parts[i].store();
      trial_move(i);
      acc_move(i, eno);
    }
  void calc_acceptance_and_adjust(void)
    {
      ntype r_tra, r_vol;
      if (tot_tra > 0)
        {
          // rate di accettazione delle mosse
          r_tra = ((double)(tot_tra - tra_rej))/tot_tra;
          std::cout << "rate tra: " << r_tra << " deltra=" << pars.deltra << "\n";
          if (r_tra > 0.5)
            {
              pars.deltra *= 1.1;
            }
          else
            {
              pars.deltra /= 1.1;
            }
          tot_tra=tra_rej = 0;
        }
      if (tot_vol > 0)
        {
          r_vol = ((ntype)(tot_vol - vol_rej)) / tot_vol;
          std::cout << "rate vol: " << r_vol << " vmax=" << pars.vmax << "\n";
          if (r_vol > 0.5) 
            {
              pars.vmax *= 1.1;
           }
          else
            {
              pars.vmax /= 1.1;
            }
          tot_vol=vol_rej=0.0;
        }
    }

  void restore_all_pars(void)
    {
      for (auto i=0; i < pars.Np; i++)
        {
          parts[i].restore();
        }
    }

  void store_all_pars(void)
    {
      for (auto i=0; i < pars.Np; i++)
        {
          parts[i].store();
        }
    }
  
  void trial_move_box(ntype& DG, ntype& fact)
    {
      ntype beta, delV, Vo, Vn, enn, eno;
      delV = pars.vmax*(rng.ranf()*2.0-1.0);
      eno = totenergy();
      Vo = pars.L[0]*pars.L[1]*pars.L[2];
      Vn = Vo + delV;
      fact = cbrt(Vn/Vo); // (Vn/Vo)^(1/3)
      for (auto i=0; i < pars.Np; i++)
        {
          parts[i].r *= fact;
        }
      pars.L *= fact;
      enn = totenergy();
      beta = 1.0/pars.T;
      DG = pars.P*delV + (enn-eno) - (1.0/beta)*pars.Np*log(Vn/Vo); 
    }
              
  void acc_move_box(ntype DG, ntype fact)
    {
      ntype beta, expDG, xi;
      beta=1.0/pars.T;
      expDG = exp(-beta*DG);
      xi = rng.ranf();
      if (xi >= expDG)
        {
          // reject 
          restore_all_pars();
          pars.L /= fact;
          vol_rej++;
        }

    }

  void move_box(void)
    {
      ntype DG, fact;
      store_all_pars();
      trial_move_box(DG, fact);
      acc_move_box(DG, fact);
    }

  void init_measures(void)
    {
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::trunc);
      f.close();
      if (pars.simtype==1)
        {
          f.open("density.dat", std::ios::out|std::ios::trunc);
          f.close();
        }
    }

 void save_measures(long int t)
    {
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::app);
      f << t << " " << totenergy()/pars.Np << "\n";
      f.close();
      if (pars.simtype==1)
        {
          f.open("density.dat", std::ios::out|std::ios::app);
          f << t << " " << pars.Np/(pars.L[0]*pars.L[1]*pars.L[2]) << "\n";
          f.close();
        }
    }
 public:
  void run(void) 
    {
      // ciclo sui passi MC
      int i, t, ip;
      tot_tra = tot_vol = tra_rej = vol_rej = 0;
      init_measures();
      for (t = 0; t < pars.totsteps; t++)
        {
          // ogni passo MC sono Np tentativi di mossa
          for (i=0; i < pars.Np; i++)
            {
              // scelgo una particella a caso
              if (pars.simtype==1 && rng.ranf() < 1.0/(pars.Np+1))
                {
                  move_box();
                  tot_vol++; 
                }
              else 
                {
                  ip = rng.ranf()*pars.Np;
                  move_NTV(ip);
                  tot_tra++;
                }
            }
          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }
          if (t > 0 && t % pars.outstps == 0) 
            {
              std::cout << "Step #" << t << "\n";
              // per confrontarsi con il Johnson si deve calcolare l'energia interna ovvero sommare il contributo
              // di energia cinetica
              std::cout << "total energy per particle is " << totenergy()/pars.Np << "\n";
              std::cout << "box size: " << pars.L[0] << " " << pars.L[1] << " " << pars.L[2] << "\n";
              if (pars.simtype==1)
                {
                  std::cout << "density is " << pars.Np / (pars.L[0]*pars.L[1]*pars.L[2]) << "\n";
                }
            }
          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }
          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && 
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }
    } 
};
template <typename ntype>

class mdsimLJ: public simLJ<ntype>
{
  using bc = simLJ<ntype>;
  using bc::pars;
  using bc::save_mgl_snapshot;
  ntype Us, U, P;

  /*Q (the "virtual mass") controls the coupling strength between the system and the thermostat. 
  A properly chosen Q ensures that the temperature fluctuates around the desired value,
  without being too erratic or too sluggish. It's a parameter that 
  requires tuning based on the specific system and conditions being simulated.*/
  ntype Q = 10; // "virtual mass" associated with the thermostat 
  /*A typical recommendation is to set  Q to a value that matches the order of magnitude 
  of the system's kinetic energy or to adjust  Q based on trial simulations to achieve stable temperature control*/

  /*xi is the thermostat variable (sometimes interpreted as the velocity of a virtual particle connected to the system). 
  It adjusts the system's kinetic energy to match the desired temperature, acting as a feedback mechanism. 
  During the simulation, xi is dynamically updated to either absorb excess kinetic energy from the system (cooling it down) 
  or add kinetic energy to the system (heating it up), thereby controlling the temperature.*/
  ntype xi = 0.0; // Initializing xi

 
  void init_measures(void)
    {
      std::fstream f;
      f.open("totenergy.dat", std::ios::out|std::ios::trunc);
      f.close();
    }
 
  void save_measures(long int t)
    {
      std::fstream f;
      f.open("totenergy.dat", std::ios::out|std::ios::app);
      ntype K=calcK();
      // K+Us is the total conserved energy, where Us is the shifted potential energy
      f << t << " " << std::setprecision(15) << K+Us << "\n";
      f.sync();      
      f.close();
    } 

  ntype calcK(void)
  {
    ntype K = 0.0;
    // [TODO] calc total kinetic energy
    for(auto &part : parts) {
        ntype v2 = part.v * part.v; // Use of the dot product for v^2 in order to calculate the square of the velocity vector's norm
        K += 0.5 * part.m * v2; // part.m is the mass of the particle
    }
    return K;
  }

  void calc_forces(ntype& Us, ntype& U)
  {
    // [TODO]
     // Us is the shifted LJ potential energy  
     // U is the cut (unshifted) LJ potential energy
     // Total force acting on each particle must be stored in parts[i].f
    Us = U = 0.0;

    // Reset forces for all particles:
    for(auto &part : parts) {
    part.f = {0, 0, 0};
    }
    
    for(int i = 0; i < parts.size() - 1; ++i) {
        for(int j = i + 1; j < parts.size(); ++j) { //grazie a questo non divido per due l'energia alla fine!
            pvector<ntype,3> fijVector;
            ntype vij, vijs, wij;
            fijVector = parts[i].fij(parts[j], pars.L, vij, vijs, wij);
            parts[i].f += fijVector; // Accumulate force on particle i
            parts[j].f -= fijVector; // Newton's third law, equal and opposite force

            U += vij; // Update unshifted LJ potential energy
            Us += vijs; // Update shifted LJ potential energy
        }
    }

  }

  void calc_pressure(ntype& P) 
  {
    ntype wij_sum = 0.0; //viriale
    ntype vij, vijs, wij;

    // Calcola la somma di wij per tutte le coppie
    for (size_t i = 0; i < parts.size(); ++i) {
      for (size_t j = i + 1; j < parts.size(); ++j) 
        {
          pvector<ntype, 3> fijv = parts[i].fij(parts[j], pars.L, vij, vijs, wij);
          wij_sum += wij;
        }
    }
    // Calcolo di P_tail
    ntype P_tail = (32.0 / 9.0) * M_PI * std::pow(pars.rho, 2) * (std::pow(pars.sigma / pars.rc, 9) - 1.5 * std::pow(pars.sigma / pars.rc, 3));
    // Calcolo della pressione totale:
    ntype volume = pars.L[0] * pars.L[1] * pars.L[2]; // Calcola il volume del sistema
    P = (pars.Np * pars.T + wij_sum / 3.0) / volume + P_tail;
    // NOTA: divido viriale per 3.0 perché - come fatto per l'energia - faccio la somma per i<j!
  }

  // Generate gaussian number by Box Muller algorithm  
  ntype gauss(void)
    {
      // gaussian of variance 1 and 0 mean
      double x1, x2;
      do
        {
          x1 = rng.ranf();
          x2 = rng.ranf();
        }
      while (x1==0 || x2==0);
      return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
    }
protected:
  using bc::parts, bc::pbc;

public:
void prepare_initial_conf(void)
  {
    bc::prepare_initial_conf(); // Initialize positions
    pvector<ntype, 3> totalVelocity;
    totalVelocity = {0, 0, 0};


    // 1) using gauss() method (Box-Muller algorithm) to initilize particle's velocity from Maxwell-Boltzmann distribution,
    for(auto &part : parts) {
        pvector<ntype, 3> vel;
        for (int d = 0; d < 3; ++d) {
            vel[d] = sqrt(pars.T / part.m) * gauss(); // Use the Box-Muller method here
        }
        part.v = vel;
        totalVelocity += vel;
    }

    // 3) set total momentum to 0
double size = static_cast<double>(parts.size());
  for (int i = 0; i < 3; ++i) {
    totalVelocity[i] /= size;
    //come dire: totalVelocity = totalVelocity / parts.size(); ma non ho fatto overloading per questo...
  }
    for(auto &part : parts) {
        part.v -= totalVelocity;
    }
    //2) set also vcut by calling parts[i].set_vcut() method
    for(auto &part : parts) {
        part.set_vcut(); 
    }
}

void run(void)
  {
    std::vector<double> energie; // Dichiarazione del vettore delle energie 
    std::vector<double> pressioni; // Dichiarazione del vettore delle pressioni 
    std::vector<double> energiePotPerParticella; // Dichiarazione del vettore dell'energiea potenziale
    calc_forces(Us, U);
    //calc_pressure(P);
    ntype energiaPotTotalePerParticella = U/pars.Np; //definisco energia potenziale per particella (da confrontarne la media U in letteratura)
    std::cout << "Initial Excess energy=" << Us << "\n";
    std::cout << "Initial Kinetic energy=" << calcK() << "\n";
  
    std::fstream fo;
    fo.open("totenergy.dat", std::fstream::out); 
    for (long long int t=0; t < pars.totsteps; t++)
      {
        // exp(iL*dt)=exp(iLp*dt/2)*exp(iLq*dt)*exp(iLp*dt/2)
        // 1) Apply exp(iLp*dt/2) to all particles
        for(auto &part : parts) {
          part.updLp(pars.dt / 2);
        }
        // 2) Applying exp(iLq*dt) to all particles and PBC:
        for (std::size_t i = 0; i < parts.size(); ++i) {
          parts[i].updLq(pars.dt); // Update position
          pbc(i); // Adjust position according to PBC
        }

        // Recalculate forces & pressure, after position update
        calc_forces(Us, U);
        energiaPotTotalePerParticella = U / pars.Np;
        // Apply VELOCITY RESCALING, only during equilibration:
        /*if (t < pars.eqstps) {
            velocityRescaling();
        }*/
        //Apply thermostat:
        NoseHooverThermostat(xi, Q);
        calc_pressure(P);

        // 3) [TODO] applying exp(iLp*dt/2) to all particles here:
        for(auto &part : parts) {
          part.updLp(pars.dt / 2);
        }

        if (pars.outstps > 0 && t % pars.outstps== 0) 
          {
            ntype K=calcK(); //calcolo energia cinetica
            std::cout << std::setprecision(6) << "[step #" << t << "] \nTotal energy = " << K+Us << "\n";
            std::cout << std::setprecision(5) << "Kinetic energy = " << K << "\n";
            std::cout << std::setprecision(5) << "Pot energy/Np = " << U/pars.Np << "\n";
            std::cout << std::setprecision(5) << "U_cs/Np = " << Us/pars.Np << "\n";
            std::cout << std::setprecision(5) << "Excess Pressure = " << P << "\n";
        }

        if (t > pars.eqstps && pars.savemeasure > 0 && (t % pars.savemeasure == 0))
          {
            save_measures(t);
            //A partire da un tempo t>tempo_di_equilibratura, eseguo i seguenti comandi:
            ntype K=calcK(); //ri-calcolo energia cinetica
            ntype energiaTotale = K + Us; // (1) Calcolo l'energia totale per il passo corrente
            energie.push_back(energiaTotale); // (2) Salvo l'energia totale nel vettore apposito
            pressioni.push_back(P); // (2) Salvo il valore della pressione nel vettore apposito
            //ntype energiaPotTotalePerParticella = U / pars.Np; //definisco energia potenziale per particella (da confrontarne la media U in letteratura)
            energiePotPerParticella.push_back(energiaPotTotalePerParticella); // memorizz valore dell'energia potenziale tot per part per il passo corrente, nel vettore
          }
      
       /* //MGL_snapshots:
       if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }*/
            if (t==9800){ 
              StampaConfig(t);
            }
         
      }
    // Calcolo e salvo la deviazione standard dell'ENERGIA TOTALE a fine simulazione
    //salvaDeviazioneStandardEnergia(energie);
    // Calcolo e salvo la deviazione standard dell'ENERGIA (potenziale) "in eccesso" PER PARTICELLA (+ tail correciton) 
    salvaEnergiaPotFinale(energiePotPerParticella);
    // Calcolo e salvo la media e la deviazione standard della pressione a fine simulazione
    salvaRisultatiPressione(pressioni);
  }

void salvaDeviazioneStandardEnergia(const std::vector<double>& energie) {
    if(energie.empty()) return; // Controlla che il vettore non sia vuoto

    double media = std::accumulate(energie.begin(), energie.end(), 0.0) / energie.size();
    double sommaQuad = std::accumulate(energie.begin(), energie.end(), 0.0, [media](double accum, double val) {
        return accum + (val - media) * (val - media);
    });
    double deviazioneStandardE = std::sqrt(sommaQuad / energie.size());
    double mediaE = media;
    // Crea il nome del file basato su pars.dt
    std::ostringstream nomeFile;
    nomeFile << std::fixed << std::setprecision(3) << "DEvsdt=" << pars.dt << ".dat";
    
    // Apri il file e scrivi la deviazione standard
    std::ofstream file(nomeFile.str());
    if(file.is_open()) {        
        file << "Deviazione Standard dell'Energia: " << deviazioneStandardE << "\n";
        file << "Passo di integrazione usato: " << pars.dt << "\n";
        file << "Step totali di integrazione: " << pars.totsteps << "\n";
        file.close();
    } else {
        std::cerr << "Impossibile aprire il file " << nomeFile.str() << std::endl;
    }
}

void salvaRisultatiPressione(const std::vector<double>& pressioni) {
    //correzione di coda alla pressione:
    if(pressioni.empty()) return;

    double media = std::accumulate(pressioni.begin(), pressioni.end(), 0.0) / pressioni.size();
    double sommaQuad = std::accumulate(pressioni.begin(), pressioni.end(), 0.0, [media](double accum, double val) {
        return accum + (val - media) * (val - media);
    });
    double deviazioneStandard = std::sqrt(sommaQuad / pressioni.size());

    // Apri il file e scrivi la media della pressione e la deviazione standard
    std::ofstream file("P*.dat");
    if(file.is_open()) {
        file << "Media Pressione + correzioni di coda: " << media << "\n";
        file << "Deviazione Standard Pressione: " << deviazioneStandard << "\n";
        file.close();
    } else {
        std::cerr << "Impossibile aprire il file P*.dat per scrivere i risultati.\n";
    }
}

void salvaEnergiaPotFinale(const std::vector<double>& energiePotPerParticella) {
    if(energiePotPerParticella.empty()) return;

    // Calcolo la correzione di coda U_tail
    double U_tail = (8.0 / 9.0) * M_PI * pars.rho * (pow(pars.sigma / pars.rc, 9) - 3 * pow(pars.sigma / pars.rc, 3));
    //Calcolo la media dell'energia potenziale per particella
    double media = std::accumulate(energiePotPerParticella.begin(), energiePotPerParticella.end(), 0.0) / energiePotPerParticella.size();
    // Aggiungo la correzione di coda alla media
    media += U_tail;

    double sommaQuad = std::accumulate(energiePotPerParticella.begin(), energiePotPerParticella.end(), 0.0, [media, U_tail](double accum, double val) {
        return accum + ((val + U_tail) - media) * ((val + U_tail) - media);
    });
    double deviazioneStandard = std::sqrt(sommaQuad / energiePotPerParticella.size());

    std::ofstream file("U*.dat");
    if(file.is_open()) {
        file << "Media Energia Potenziale per Particella (con correzione di coda): " << media << "\n";
        file << "Deviazione Standard Energia Potenziale per Particella: " << deviazioneStandard << "\n";
        file.close();
    } else {
        std::cerr << "Impossibile aprire il file U*.dat per scrivere i risultati.\n";
    }
}

void velocityRescaling() {
  double currentTemperature = calcCurrentTemperature(); 
  double targetTemperature = pars.T; // The target temperature of my system
  double rescalingFactor = sqrt(targetTemperature / currentTemperature);

  for (auto &part : parts) {
    part.v *= rescalingFactor; // Rescale velocities
  }
}

ntype calcCurrentTemperature() {
  ntype kineticEnergy = 0.0;
  for (const auto &particle : parts) {
    kineticEnergy += 0.5 * particle.m * (particle.v * particle.v); 
  }

  ntype degreesOfFreedom = 3 * parts.size();
  ntype currentTemperature = (2.0 / degreesOfFreedom) * kineticEnergy; // Equipartition theorem: (2/2N)*k_B*K

  return currentTemperature;
}

// Function to update velocities and xi for Nosé-Hoover thermostat
void NoseHooverThermostat(ntype& xi, ntype Q) {
    ntype kineticEnergy = calcK();
    ntype targetTemperature = pars.T; // The target temperature of my system
    ntype gkT = (3 * parts.size() + 1) * targetTemperature; // gkT is the target kinetic energy, with g=3N+1 and k_b=1
    xi += pars.dt * ((2.0 * kineticEnergy - gkT) / Q); // Update xi
    
    // Update velocities based on xi
  for (auto &part : parts) {
    part.v += -xi * pars.dt * part.v;
    }
}


void StampaConfig(long long int t)
{
    std::fstream f;
    std::string s;

    s = "config-" + std::to_string(t) + ".dat";

    f.open(s, std::ios::out | std::ios::trunc);
    for (int i = 0; i < pars.Np; i++)
    {
        f << parts[i].r[0] << " " << parts[i].r[1] << " " << parts[i].r[2] << "\n";
        //f << parts[i].v[0] << " " << parts[i].v[1] << " " << parts[i].v[2] << "\n";
    }
    f.close();
}



};
#endif
