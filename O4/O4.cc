#include "O4.h"

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
  int MPO_cutoff = 1E-20;
  
  auto N = input.getInt("N");
  bool periodic = input.getInt("periodic");
  auto a1 = input.getReal("a1");
  auto a2 = input.getReal("a2");
  auto a3 = input.getReal("a3");
  auto a4 = input.getReal("a4");
  auto U = input.getReal("U");
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
  bool measure_dimerVEV = input.getInt("measure_dimerVEV");
  bool measure_dimer2 = input.getInt("measure_dimer2");
  
  /* Reuse previous mps */
  Electron sites_ini = Electron(N);
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
      auto a2_ini = input.getReal("a2_ini");
      auto a3_ini = input.getReal("a3_ini");
      auto U_ini = input.getReal("U_ini");
      auto dmrg_cutoff_ini = input.getReal("dmrg_cutoff_ini");
      
      auto sites_ini_file = format("data_mps/sites_periodic%d_size%d_a2%.1f_a3%.1f_U%.5f_dmrg_cutoff%.1E.dat", periodic, N, a2_ini, a3_ini, U_ini, dmrg_cutoff_ini);
      readFromFile(sites_ini_file, sites_ini);
      
      for (int i = 0; i < nstates; ++i)
	{
	  auto psi_ini_file = format("data_mps/mps_periodic%d_size%d_a2%.1f_a3%.1f_U%.5f_dmrg_cutoff%.1E_state%d.dat", periodic, N, a2_ini, a3_ini, U_ini, dmrg_cutoff_ini, i);
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
  
  O4 O4(N, periodic, a1, a2, a3, a4, U, MPO_cutoff);
  O4.sites = sites_ini;
  O4.H = O4.O4Hamiltonian();

  /*auto sw_group = InputGroup(argv[1],"sweep_group");
  auto sw_table = InputGroup(sw_group,"sweep_table");
  auto sw_table1 = InputGroup(sw_group,"sweep_table1");
  auto ini_sweeps = Sweeps(ini_nsweep,sw_table1);
  auto sweeps = Sweeps(nsweep,sw_table1);*/

  Measurements GroundStates(new_measurements, dist_to_boundary, precision, ini_nsweep, nsweep, ini_sweeps, sweeps, penalty_weight);
  std::vector<Real> en(nstates);
  std::vector<MPS> psi(nstates);

  GroundStates.sites_file = format("data_mps/sites_periodic%d_size%d_a2%.1f_a3%.1f_U%.5f_dmrg_cutoff%.1E.dat", periodic, N, O4.a2, O4.a3, O4.U, dmrg_cutoff);
  GroundStates.energy_file = format("data_mps/energy_periodic%d_size%d_a2%.1f_a3%.1f_U%.5f_dmrg_cutoff%.1E.dat", periodic, N, O4.a2, O4.a3, O4.U, dmrg_cutoff);
  GroundStates.mps_file.resize(nstates);
  for (int i = 0; i < nstates; ++i)
    {
      GroundStates.mps_file[i] = format("data_mps/mps_periodic%d_size%d_a2%.1f_a3%.1f_U%.5f_dmrg_cutoff%.1E_state%d.dat", periodic, O4.N, O4.a2, O4.a3, O4.U, dmrg_cutoff, i);
    }

  if (!utils::fileExists(GroundStates.mps_file[nstates-1]) || new_dmrg)
    {
      writeToFile(GroundStates.sites_file, O4.sites);
      GroundStates.GroundStates(O4, nstates, psi_ini, en, psi);
    }
  else
    {
      readFromFile(GroundStates.sites_file, O4.sites);
      for (int i = 0; i < nstates; ++i)
	readFromFile(GroundStates.mps_file[i], psi[i]);
    }

  printfln("N = %d\t periodic = %d\t a2 = %.3f\t a3 = %.3f\t U = %.5f dmrg cutoff = %.1E\t precision = %.1E\t penalty weight = %d", N, periodic, a2, a3, U, dmrg_cutoff, precision, penalty_weight);
  
  if (measure_spin) GroundStates.PrintSpin(O4, psi[0]);
  if (measure_dimer) GroundStates.PrintDimer(O4, psi[0]);
  if (measure_dimerVEV) GroundStates.PrintDimerVEV(O4, psi[0]);
  //GroundStates.PrintfSpinxCorrelation(O4, psi[0]);
  //GroundStates.PrintfEntanglementEntropy(O4, psi[0]);

  //GroundStates.PrintMomentum(O4, psi[0]);
  
  return 0;
}

