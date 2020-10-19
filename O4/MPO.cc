#include "O4.h"

using namespace itensor;

MPO Measurements::S2(Electron const& sites)
{
  auto N = length(sites);
  
  auto ampo = AutoMPO(sites);

  for(int i = 1; i <= N; ++i)
    {
      for(int j = 1; j <= N; ++j)
	{
	  ampo += 0.5,"S+",i,"S-",j;
	  ampo += 0.5,"S-",i,"S+",j;
	  ampo += 1.0,"Sz",i,"Sz",j;
	}
    }
  
  auto S2 = toMPO(ampo,{"Cutoff=",1E-20});
  
  return S2;
}



MPO Measurements::Momentum(Electron const& sites)
{
  auto N = length(sites);
  
  auto ampo = AutoMPO(sites);
  //std::complex<double> I(0.0, 1.0);
  
  for(int k = 1; k <= N; ++k)
    {
      for(int i = 1; i <= N; ++i)
	{
	  for(int j = 1; j <= N; ++j)
	    {
	      ampo += 2.0*PI/(N*N)*k*std::exp(Cplx_i*2.0*PI/N*k*(i-j)),"Cdagup",j,"Cup",i;
	      ampo += 2.0*PI/(N*N)*k*std::exp(Cplx_i*2.0*PI/N*k*(i-j)),"Cdagdn",j,"Cdn",i;
	    }
	}
    }

  //ampo += 1;
  
  auto P = toMPO(ampo,{"Cutoff=",1E-20});
  //auto T = toExpH(ampo,0.1,{"Cutoff=",1E-20});
  
  return P;
}
