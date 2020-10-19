#ifndef SpinChain_H
#define SpinChain_H

#include "itensor/all.h"

using namespace itensor;

class Hamiltonian
{
public:
  int N;
  bool periodic;
  double J1;
  double J2;
  double cutoff;
  SpinHalf sites;
  MPO H;
  
  Hamiltonian(int N, bool periodic, double J1, double J2, double cutoff): N(N), periodic(periodic), J1(J1), J2(J2), cutoff(cutoff) {};
  
  MPO SpinChain();
};

class Measurements
{
public:
  bool new_measurements;
  int dist_to_boundary;
  double precision;
  int ini_nsweep;
  int nsweep;
  Sweeps ini_sweeps;
  Sweeps sweeps;
  double penalty_weight;
  string sites_file;
  string energy_file;
  std::vector<string> mps_file;

  std::chrono::high_resolution_clock::time_point start;
  
  Measurements(bool new_measurements, int dist_to_boundary, double precision, int ini_nsweep, int nsweep, Sweeps ini_sweeps, Sweeps sweeps, double penalty_weight): new_measurements(new_measurements), dist_to_boundary(dist_to_boundary), precision(precision), ini_nsweep(ini_nsweep), nsweep(nsweep), ini_sweeps(ini_sweeps), sweeps(sweeps), penalty_weight(penalty_weight), start(std::chrono::high_resolution_clock::now()) {};

  void GroundStates(const Hamiltonian&, int, const std::vector<MPS>&, std::vector<Real>&, std::vector<MPS>&);
  Real TwopointCorrelation(const SpinHalf&, MPS&, string, string, int, int);
  Real FourpointCorrelation(const SpinHalf&, MPS&, string, int, int, int, int);
  Real DimerCorrelation(const SpinHalf&, MPS&, string, int, int);
  Real EntanglementEntropy(MPS&, int);
  void PrintfSpinCorrelation(const Hamiltonian&, MPS&);
  void PrintfSpinxCorrelation(const Hamiltonian&, MPS&);
  void PrintfDimerCorrelation(const Hamiltonian&, MPS&);
  void PrintfDimer2Correlation(const Hamiltonian&, MPS&);
  void PrintfEntanglementEntropy(const Hamiltonian&, MPS&);
};
  

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  if ( !v.empty() )
    {
      for (auto it = v.begin(); it != v.end(); ++it)
	out << *it << " ";
    }
  
  return out;
}

namespace utils {
  
  inline bool fileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
  }

}

#endif
