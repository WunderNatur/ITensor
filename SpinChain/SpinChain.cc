#include "SpinChain.h"

using namespace itensor;

int main(int argc, char* argv[]) 
{
  if(argc != 2) 
    { 
      //reminds us to give an input file if we forget
      printfln("Usage: %s inputfile",argv[0]); 
      return 0; 
    }
  auto input = InputGroup(argv[1],"input");

  bool new_dmrg = input.getInt("new_dmrg");
  auto MPO_cutoff = 1E-20;
  
  auto N = input.getInt("N");
  bool periodic = input.getInt("periodic");
  auto J1 = input.getReal("J1");
  auto J2 = input.getReal("J2");
  auto dmrg_cutoff = input.getReal("dmrg_cutoff");
  auto precision = input.getReal("precision");
  auto penalty_weight = input.getReal("penalty_weight");
  auto nstates = input.getInt("nstates");
  int ini_nsweep = 1;
  int nsweep = 1;
  auto sweeps = Sweeps(nsweep, 20, 1E4, dmrg_cutoff, 0.);
  auto ini_sweeps = sweeps;
  
  bool new_measurements = input.getInt("new_measurements");
  auto dist_to_boundary = input.getInt("dist_to_boundary");
  bool measure_spin = input.getInt("measure_spin");
  bool measure_dimer = input.getInt("measure_dimer");

  /* Reuse previous mps */
  SpinHalf sites_ini = SpinHalf(N);
  std::vector<MPS> psi_ini(nstates);

  auto state = InitState(sites_ini);
  for(auto i : range1(N))
    {
      if(i%2 == 1) state.set(i,"Up");
      else         state.set(i,"Dn");
    }
  auto psi_ran = randomMPS(state);

  bool use_psi_ini = input.getInt("use_psi_ini");
  if (use_psi_ini)
    {
      auto J2_ini = input.getReal("J2_ini");
      auto dmrg_cutoff_ini = input.getReal("dmrg_cutoff_ini");
      
      auto sites_ini_file = format("data_mps/sites_periodic%d_size%d_J2%.5f_dmrg_cutoff%.1E.dat", periodic, N, J2_ini, dmrg_cutoff_ini);
      readFromFile(sites_ini_file, sites_ini);
      
      for (int i = 0; i < nstates; ++i)
	{
	  auto psi_ini_file = format("data_mps/mps_periodic%d_size%d_J2%.5f_dmrg_cutoff%.1E_state%d.dat", periodic, N, J2_ini, dmrg_cutoff_ini, i);
	  if (utils::fileExists(psi_ini_file))
	    readFromFile(psi_ini_file, psi_ini[i]);
	  else
	    psi_ini[i] = psi_ran;
	}
    }
  else
    {
      for (int i = 0; i < nstates; ++i) psi_ini[i] = psi_ran;
    }

  Hamiltonian H(N, periodic, J1, J2, MPO_cutoff);
  H.sites = sites_ini;
  H.H = H.SpinChain();

  /*auto sw_group = InputGroup(argv[1],"sweep_group");
  auto sw_table = InputGroup(sw_group,"sweep_table");
  auto sw_table1 = InputGroup(sw_group,"sweep_table1");
  auto ini_sweeps = Sweeps(ini_nsweep,sw_table1);
  auto sweeps = Sweeps(nsweep,sw_table1);*/

  Measurements GroundStates(new_measurements, dist_to_boundary, precision, ini_nsweep, nsweep, ini_sweeps, sweeps, penalty_weight);
  std::vector<Real> en(nstates);
  std::vector<MPS> psi(nstates);

  GroundStates.sites_file = format("data_mps/sites_periodic%d_size%d_J2%.5f_dmrg_cutoff%.1E.dat", periodic, N, J2, dmrg_cutoff);
  GroundStates.energy_file = format("data_mps/energy_periodic%d_size%d_J2%.5f_dmrg_cutoff%.1E.dat", periodic, N, J2, dmrg_cutoff);
  GroundStates.mps_file.resize(nstates);
  for (int i = 0; i < nstates; ++i)
    {
      GroundStates.mps_file[i] = format("data_mps/mps_periodic%d_size%d_J2%.5f_dmrg_cutoff%.1E_state%d.dat", periodic, N, J2, dmrg_cutoff, i);
    }

  if (!utils::fileExists(GroundStates.mps_file[nstates-1]) || new_dmrg)
    {
      writeToFile(GroundStates.sites_file, H.sites);
      GroundStates.GroundStates(H, nstates, psi_ini, en, psi);
    }
  else
    {
      readFromFile(GroundStates.sites_file, H.sites);
      for (int i = 0; i < nstates; ++i)
	readFromFile(GroundStates.mps_file[i], psi[i]);
    }

  printfln("N = %d\t periodic = %d\t J2 = %.5f\t dmrg cutoff = %.1E\t precision = %.1E\t penalty weight = %d", N, periodic, J2, dmrg_cutoff, precision, penalty_weight);
  
  if (measure_spin) GroundStates.PrintfSpinCorrelation(H, psi[0]);
  if (measure_dimer) GroundStates.PrintfDimerCorrelation(H, psi[0]);
  //GroundStates.PrintfSpinxCorrelation(H, psi[0]);
  //GroundStates.PrintfEntanglementEntropy(H, psi[0]);

  //std::cout << H.H;
  
  return 0;
}

