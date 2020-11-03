#include "O4.h"

using namespace itensor;

void Measurements::GroundStates(const O4& O4, int n, const std::vector<MPS>& psi_ini, std::vector<Real>& en, std::vector<MPS>& psi)
{
  int sweep_count_total = 0;
  std::vector<int> sweep_count(n,0);
  auto wfs = std::vector<MPS>(0);
  
  std::ofstream pfile;
  pfile.open (energy_file);
  pfile.precision(15);
      
  for (int i = 0; i < n; ++i)
    {
      auto [en_,psi_] = dmrg(O4.H, wfs, psi_ini[i], ini_sweeps,{"Quiet",true,"Weight=",penalty_weight});
      en[i] = en_;
      psi[i] = psi_;
      sweep_count[i] += ini_nsweep;
      sweep_count_total += ini_nsweep;

      double en_change = 1;
      while(abs(en_change) > precision)
	{
	  auto [en_,psi_] = dmrg(O4.H, wfs, psi[i],sweeps,{"Quiet",true,"Weight=",penalty_weight});
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

Real Measurements::TwopointCorrelation(const Electron& sites, MPS& psi, string Opi, string Opj, int site_i, int site_j)
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

Real Measurements::FourpointCorrelation(const Electron& sites, MPS& psi, string Op, int site_i, int site_j, int site_k, int site_l)
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

Real Measurements::DimerCorrelation(const Electron& sites, MPS& psi, string Op, int site_i, int site_j)
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



void Measurements::PrintfSpinCorrelation(const O4& O4, MPS& psi)
{
  auto N = O4.N;
  auto file = format("data/spin_size%d_a2%.1f_a3%.1f_U%.5f.dat", N, O4.a2, O4.a3, O4.U);
  //println(psi);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = O4.sites;
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

void Measurements::PrintfSpinxCorrelation(const O4& O4, MPS& psi)
{
  auto N = O4.N;
  auto file = format("data/spinx_size%d_a2%.1f_a3%.1f_U%.5f.dat", N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = O4.sites;
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

void Measurements::PrintfDimerCorrelation(const O4& O4, MPS& psi)
{
  auto N = O4.N;
  auto file = format("data/dimer_size%d_a2%.1f_a3%.1f_U%.5f.dat", N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = O4.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      int max_dist = N/2-dist_to_boundary-1;
      std::vector<Real> corN(max_dist);
      std::vector<Real> corN_1(max_dist);
      for (int i = 0; i <= max_dist - 1; ++i)
	{
	  auto corN_1i = DimerCorrelation(sites, psi, "Sz", N/2-2, N/2+i+1);
	  corN_1[i] = corN_1i;
	  auto corNi = DimerCorrelation(sites, psi, "Sz", N/2-1, N/2+i+1);
	  corN[i] = corNi;
	}
      
      for (int i = 3; i <= max_dist+1; ++i)
	{
	  ofile << i << "\t" << std::pow(1, i) * (corN_1[i-3]-corN_1[i-2]-corN[i-3]+corN[i-2])/4.0 << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Dimer time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfDimer2Correlation(const O4& O4, MPS& psi)
{
  auto N = O4.N;
  auto file = format("data/dimer2_size%d_a2%.1f_a3%.1f_U%.5f.dat", N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      auto sites = O4.sites;
      std::ofstream ofile;
      ofile.open (file);
      
      int max_dist = N/2-dist_to_boundary-1;
      for (int i = 2; i<= max_dist; ++i)
	{
	  auto cor = FourpointCorrelation(sites, psi, "Sz", N/2, N/2+1, N/2+i, N/2+i+1)*16;
	  ofile << i << "\t" << cor << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Dimer2 time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintfEntanglementEntropy(const O4& O4, MPS& psi)
{
  auto file = format("data/entropy_size%d_a2%.1f_a3%.1f_U%.5f.dat", O4.N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      std::ofstream ofile;
      ofile.open (file);
      
      for (int i = 1; i<O4.N; ++i)
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

void Measurements::PrintS2(const O4& O4, MPS& psi)
{
  auto S2mpo = S2(O4.sites);
  auto S2num = inner(psi, S2mpo, psi);
  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  printfln("Total spin = %.5f", S2num);
  std::cout << "Total spin time = " << duration.count() << " seconds\n";
}

void Measurements::PrintMomentum(const O4& O4, MPS& psi)
{
  auto Pmpo1 = Momentum(O4.sites);
  auto Pmpo2 = nmultMPO(Pmpo1, prime(Pmpo1));
  auto Pmpo3 = nmultMPO(Pmpo1, prime(Pmpo2));
  auto Pmpo4 = nmultMPO(Pmpo1, prime(Pmpo3));
  auto Pmpo5 = nmultMPO(Pmpo1, prime(Pmpo4));
  auto Pmpo6 = nmultMPO(Pmpo1, prime(Pmpo5));
  auto Pnum1 = innerC(psi, Pmpo1, psi);
  auto Pnum2 = innerC(psi, Pmpo2, psi);
  auto Pnum3 = innerC(psi, Pmpo3, psi);
  auto Pnum4 = innerC(psi, Pmpo4, psi);
  auto Pnum5 = innerC(psi, Pmpo5, psi);
  auto Pnum6 = innerC(psi, Pmpo6, psi);
  auto Pnum = 1 + Pnum1 + Pnum2/2.0 + Pnum3/6.0 + Pnum4/24.0 + Pnum5/120.0 + Pnum6/720.0;
  auto Pnum_r = fmod(real(Pnum),2);
  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  printfln("P = %.5f pi", Pnum);
  printfln("Momentum = %.5f pi", Pnum_r);
  std::cout << "Momentum time = " << duration.count() << " seconds\n";
}

void Measurements::PrintSpin(const O4& O4, MPS& psi)
{
  auto site_i = dist_to_boundary;
  auto file = format("data/spin_size%d_a2%.1f_a3%.1f_U%.5f.dat", O4.N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      std::ofstream ofile;
      ofile.open (file);
      for (int i=0; i<=O4.N-site_i; ++i)
	{
	  auto spin_mpo = Spin(O4.sites, site_i, site_i + i);
	  auto spin_num = inner(psi, spin_mpo, psi);
	  
	  ofile << i << "\t" << std::pow(-1, i) * spin_num << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Spin time = " << duration.count() << " seconds\n";
    }
}

void Measurements::PrintDimer(const O4& O4, MPS& psi)
{
  auto site_i = dist_to_boundary;
  auto file = format("data/dimer_size%d_a2%.1f_a3%.1f_U%.5f.dat", O4.N, O4.a2, O4.a3, O4.U);
  if (!utils::fileExists(file) || new_measurements)
    {
      std::ofstream ofile;
      ofile.open (file);
      for (int i=0; i<=O4.N - site_i -1; ++i)
	{
	  auto dimer_mpo = Dimer(O4.sites, site_i, site_i+1, site_i + i, (site_i + i)%O4.N+1);
	  auto dimer_num = inner(psi, dimer_mpo, psi);
	  
	  ofile << i << "\t" << std::pow(1, i) * dimer_num << "\n";
	}
      
      ofile.close();
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Dimer time = " << duration.count() << " seconds\n";
    }
}
