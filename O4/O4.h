#ifndef O4_H
#define O4_H

#include "itensor/all.h"

using namespace itensor;

# define PI 3.141592653589793238462643383279502884

class O4
{
public:
  int N;
  bool periodic;
  double a1;
  double a2;
  double a3;
  double a4;
  double U;
  double cutoff;
  Electron sites;
  //InitState state;
  MPO H;
  
  O4(int N, bool periodic, double a1, double a2, double a3, double a4, double U, double cutoff): N(N), periodic(periodic), a1(a1), a2(a2), a3(a3), a4(a4), U(U), cutoff(cutoff) {};
  
  MPO O4Hamiltonian();
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

  void GroundStates(const O4&, int, const std::vector<MPS>&, std::vector<Real>&, std::vector<MPS>&);
  Real TwopointCorrelation(const Electron&, MPS&, string, string, int, int);
  Real FourpointCorrelation(const Electron&, MPS&, string, int, int, int, int);
  Real DimerCorrelation(const Electron&, MPS&, string, int, int);
  Real EntanglementEntropy(MPS&, int);
  void PrintfSpinCorrelation(const O4&, MPS&);
  void PrintfSpinxCorrelation(const O4&, MPS&);
  void PrintfDimerCorrelation(const O4&, MPS&);
  void PrintfDimer2Correlation(const O4&, MPS&);
  void PrintfEntanglementEntropy(const O4&, MPS&);

  MPO S2(Electron const&);
  MPO Momentum(Electron const&);
  MPO Spin(Electron const&, int, int);
  MPO Dimer(Electron const&, int, int, int, int);
  void PrintS2(const O4&, MPS&);
  void PrintMomentum(const O4&, MPS&);
  void PrintSpin(const O4&, MPS&);
  void PrintDimer(const O4&, MPS&);
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
