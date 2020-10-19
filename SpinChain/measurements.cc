#include "SpinChain.h"

using namespace itensor;

void Measurements::GroundStates(const Hamiltonian& H, int n, const std::vector<MPS>& psi_ini, std::vector<Real>& en, std::vector<MPS>& psi)
{
  int sweep_count_total = 0;
  std::vector<int> sweep_count(n,0);
  auto wfs = std::vector<MPS>(0);
  
  std::ofstream pfile;
  pfile.open(energy_file);
  pfile.precision(15);
  
  for (int i = 0; i < n; ++i)
    {
      auto [en_,psi_] = dmrg(H.H, wfs, psi_ini[i],ini_sweeps,{"Quiet",true,"Weight=",penalty_weight});
      en[i] = en_;
      psi[i] = psi_;
      sweep_count[i] += ini_nsweep;
      sweep_count_total += ini_nsweep;

      double en_change = 1;
      while(abs(en_change) > precision)
	{
	  auto [en_,psi_] = dmrg(H.H, wfs, psi[i],sweeps,{"Quiet",true,"Weight=",penalty_weight});
	  psi[i] = psi_;
	  en_change = en_ - en[i];
	  en[i] = en_;
	  sweep_count[i] += nsweep;
	  sweep_count_total += nsweep;
	  auto stop = std::chrono::high_resolution_clock::now();
	  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	  printfln("Ground State Energy Change = %.2E", en_change);
	  std::cout << "Sweep count = " << sweep_count << "Time taken = " << duration.count() << " seconds\n";
	}

      wfs.push_back(psi[i]);
      writeToFile(mps_file[i], psi[i]);
      pfile << en[i] << "\n";
    }

  pfile.close();
  
  for (int i = 0; i < n; ++i)
    printfln("\nen[%d] = %.10f",i,en[i]);

  for (int i = 0; i < n; ++i)
    {
      for (int j = i+1; j < n; ++j)
	printfln("\nOverlap <psi%d|psi%d> = %.2E",i,j,inner(psi[i],psi[j]));
    }
}

Real Measurements::TwopointCorrelation(const SpinHalf& sites, MPS& psi, string Opi, string Opj, int site_i, int site_j)
{
  int N = length(sites);
  
  //Make the operators you want to measure
  auto op_i = op(sites,Opi,site_i);
  auto op_j = op(sites,Opj,site_j);
  
  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(site_i); 
  
  //Create the bra/dual version of the MPS psi
  auto psidag = dag(psi);
  
  //Prime the link indices to make them distinct from
  //the original ket links
  psidag.prime("Link");
  
  //index linking i-1 to i:
  auto li_1 = leftLinkIndex(psi,site_i);
  
  auto C = prime(psi(site_i),li_1) * op_i * prime(psidag(site_i),"Site");
  
  for(int k = site_i+1; k < site_j; ++k)
    {
      C *= psi(k) * psidag(k);
    }
  //index linking j to j+1:
  auto lj = rightLinkIndex(psi,site_j);
  
  C *= prime(psi(site_j),lj) * op_j * prime(psidag(site_j),"Site");
  
  auto result = elt(C); //or eltC(C) if expecting complex

  //printfln("\nCorrelation = %.10f",4*result);
  return result;
}

Real Measurements::FourpointCorrelation(const SpinHalf& sites, MPS& psi, string Op, int site_i, int site_j, int site_k, int site_l)
{
  int N = length(sites);
  
  //Make the operators you want to measure
  auto op_i = op(sites,Op,site_i);
  auto op_j = op(sites,Op,site_j);
  auto op_k = op(sites,Op,site_k);
  auto op_l = op(sites,Op,site_l);
  
  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(site_i); 
  
  //Create the bra/dual version of the MPS psi
  auto psidag = dag(psi);
  
  //Prime the link indices to make them distinct from
  //the original ket links
  psidag.prime("Link");
  
  //index linking i-1 to i:
  auto li_1 = leftLinkIndex(psi,site_i);
  
  auto C = prime(psi(site_i),li_1) * op_i * prime(psidag(site_i),"Site");
  for(int i = site_i+1; i < site_j; ++i)
    {
      C *= psi(i) * psidag(i);
    }
  C *= psi(site_j) * op_j * prime(psidag(site_j),"Site");
  for(int i = site_j+1; i < site_k; ++i)
    {
      C *= psi(i) * psidag(i);
    }
  C *= psi(site_k) * op_k * prime(psidag(site_k),"Site");
  for(int i = site_k+1; i < site_l; ++i)
    {
      C *= psi(i) * psidag(i);
    }
  //index linking j to j+1:
  auto ll = rightLinkIndex(psi,site_l);
  
  C *= prime(psi(site_l),ll) * op_l * prime(psidag(site_l),"Site");
  
  auto result = elt(C); //or eltC(C) if expecting complex

  //printfln("\nCorrelation = %.10f",4*result);
  return result;
}

Real Measurements::DimerCorrelation(const SpinHalf& sites, MPS& psi, string Op, int site_i, int site_j)
{
  int N = length(sites);
  
  //Make the operators you want to measure
  auto op_i = op(sites,Op,site_i);
  auto op_i1 = op(sites,Op,site_i+1);
  auto op_j = op(sites,Op,site_j);
  auto op_j1 = op(sites,Op,site_j+1);
  
  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(site_i); 
  
  //Create the bra/dual version of the MPS psi
  auto psidag = dag(psi);
  
  //Prime the link indices to make them distinct from
  //the original ket links
  psidag.prime("Link");
  
  //index linking i-1 to i:
  auto li_1 = leftLinkIndex(psi,site_i);
  
  auto C = prime(psi(site_i),li_1) * op_i *  prime(psidag(site_i),"Site");
  C *= psi(site_i+1) * op_i1 * prime(psidag(site_i+1),"Site");
  for(int k = site_i+2; k < site_j; ++k)
    {
      C *= psi(k) * psidag(k);
    }
  //index linking j to j+1:
  auto lj = rightLinkIndex(psi,site_j+1);
  
  C *= psi(site_j) * op_j * prime(psidag(site_j),"Site");
  C *= prime(psi(site_j+1),lj) * op_j1 * prime(psidag(site_j+1),"Site");
  
  auto result = elt(C); //or eltC(C) if expecting complex

  //printfln("\nCorrelation = %.10f",4*result);
  return result;
}

Real Measurements::EntanglementEntropy(MPS& psi, int b)
{
  int N = length(psi);
  psi.position(b); 
  
  //SVD this wavefunction to get the spectrum
  //of density-matrix eigenvalues
  auto l = leftLinkIndex(psi,b);
  auto s = siteIndex(psi,b);
  auto [U,S,V] = svd(psi(b),{l,s});
  auto u = commonIndex(U,S);
  
  //Apply von Neumann formula
  //to the squares of the singular values
  Real SvN = 0.;
  for(auto n : range1(dim(u)))
    {
      auto Sn = elt(S,n,n);
      auto p = sqr(Sn);
      if(p > 1E-12) SvN += -p*log(p);
    }
  //printfln("Across bond b=%d, SvN = %.10f",b,SvN);

  return SvN;
}



void Measurements::PrintfSpinCorrelation(const Hamiltonian& H, MPS& psi)
{
  auto N = H.N;
  auto file = format("data/spin_size%d_J2%.5f.dat", N, H.J2);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = H.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      for (int i = 1; i<= N/2-dist_to_boundary; ++i)
	{
	  auto cor = TwopointCorrelation(sites, psi, "Sz", "Sz", N/2, N/2 + i);
	  ofile << i << "\t" << std::pow(1, i) * cor << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Spin time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfSpinxCorrelation(const Hamiltonian& H, MPS& psi)
{
  auto N = H.N;
  auto file = format("data/spinx_size%d_J2%.5f.dat", N, H.J2);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = H.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      for (int i = 1; i<= N/2-dist_to_boundary; ++i)
	{
	  auto cor = TwopointCorrelation(sites, psi, "S+", "S-", N/2, N/2 + i);
	  ofile << i << "\t" << std::pow(1, i) * cor << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Spinx time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfDimerCorrelation(const Hamiltonian& H, MPS& psi)
{
  auto N = H.N;
  auto file = format("data/dimer_size%d_J2%.5f.dat", N, H.J2);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = H.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      int max_dist = N/2-dist_to_boundary-1;
      std::vector<Real> corN(max_dist - 1);
      std::vector<Real> corN_1(max_dist - 1);
      for (int i = 0; i <= max_dist - 2; ++i)
	{
	  auto corN_1i = DimerCorrelation(sites, psi, "Sz", N/2-1, N/2+i+2);
	  corN_1[i] = corN_1i;
	  auto corNi = DimerCorrelation(sites, psi, "Sz", N/2, N/2+i+2);
	  corN[i] = corNi;
	}
      
      for (int i = 3; i <= max_dist; ++i)
	{
	  ofile << i << "\t" << std::pow(1, i) * (corN_1[i-3]-corN_1[i-2]-corN[i-3]+corN[i-2])/4.0 << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Dimer time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfDimer2Correlation(const Hamiltonian& H, MPS& psi)
{
  auto N = H.N;
  auto file = format("data/dimer2_size%d_J2%.1f.dat", N, H.J2);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = H.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      int max_dist = N/2-dist_to_boundary-2;
      for (int i = 3; i<= max_dist; ++i)
	{
	  auto cor = FourpointCorrelation(sites, psi, "Sz", N/2, N/2+2, N/2+i, N/2+i+2);
	  ofile << i << "\t" << cor << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Dimer2 time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfEntanglementEntropy(const Hamiltonian& H, MPS& psi)
{
  auto file = format("data/entropy_size%d_J2%.5f.dat", H.N, H.J2);
  if (!utils::fileExists(file) || new_measurements)
    {
      std::ofstream ofile;
      ofile.open (file);
      
      for (int i = 1; i< H.N; ++i)
	{
	  auto cor = EntanglementEntropy(psi, i);
	  ofile << i << "\t" << cor << "\n";
    }
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Entropy time = " << duration.count() << " seconds\n";
    }
}
