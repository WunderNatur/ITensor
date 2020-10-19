#include "O4_2d.h"

using namespace itensor;

MPO O4::O4Hamiltonian()
{
  auto ampoO4 = AutoMPO(sites);
  auto lattice = squareLattice(Nx,Ny,{"YPeriodic=",yperiodic});
  for(int i = 1; i <= N; ++i)
    {
      ampoO4 += U,"Nupdn",i;
      ampoO4 += -U/2,"Nup",i;
      ampoO4 += -U/2,"Ndn",i;
      ampoO4 += U/4,"Id",i;
    }

  int eta;
  for(auto bnd : lattice)
    {
      if (bnd.s1%2 && bnd.s2%2 && pi_flux) eta = -1;
      else eta = 1;
      
      ampoO4 = BondHamiltonian(ampoO4, eta, bnd.s1, bnd.s2);

      if (Ny==2 && (bnd.s1+bnd.s2)%2)
      ampoO4 = BondHamiltonian(ampoO4, eta, bnd.s1, bnd.s2);
    }

  if (xperiodic)// && Nx!=2)
    {
      for (int i = 1; i <= Ny; ++i)
	{
	  if (i%2 && (Ny*(Nx-1)+i)%2 && pi_flux) eta = -1;
	  else eta = 1;
	  
	  ampoO4 = BondHamiltonian(ampoO4, eta, Ny*(Nx-1)+i, i);
	}
    }
  
  auto H = toMPO(ampoO4,{"Cutoff=",cutoff});

  return H;
}

AutoMPO O4::BondHamiltonian(AutoMPO& ampo, int eta, int bnd_i, int bnd_j)
{
  ampo += -0.25,"Id",bnd_i;
  ampo += -eta*a1/2,"Cdagup",bnd_i,"Cup",bnd_j;
  ampo += -eta*a1/2,"Cdagup",bnd_j,"Cup",bnd_i;
  ampo += -eta*a1/2,"Cdagdn",bnd_i,"Cdn",bnd_j;
  ampo += -eta*a1/2,"Cdagdn",bnd_j,"Cdn",bnd_i;
  ampo += -eta*eta*a2,"Cdagup",bnd_i,"Cup",bnd_j,"Cdagdn",bnd_i,"Cdn",bnd_j;
  ampo += -eta*eta*a2,"Cdagup",bnd_i,"Cup",bnd_j,"Cdagdn",bnd_j,"Cdn",bnd_i;
  ampo += -eta*eta*a2,"Cdagup",bnd_j,"Cup",bnd_i,"Cdagdn",bnd_i,"Cdn",bnd_j;
  ampo += -eta*eta*a2,"Cdagup",bnd_j,"Cup",bnd_i,"Cdagdn",bnd_j,"Cdn",bnd_i;
  ampo += eta*eta*a2,"Nup",bnd_i,"Nup",bnd_j;
  ampo += -eta*eta*a2/2,"Nup",bnd_i;
  ampo += -eta*eta*a2/2,"Nup",bnd_j;
  ampo += eta*eta*a2/4,"Id",bnd_i;
  ampo += eta*eta*a2,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += -eta*eta*a2/2,"Ndn",bnd_i;
  ampo += -eta*eta*a2/2,"Ndn",bnd_j;
  ampo += eta*eta*a2/4,"Id",bnd_i;
  ampo += eta*eta*eta*2*a3,"Cdagup",bnd_i,"Cup",bnd_j,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += -eta*eta*eta*a3,"Cdagup",bnd_i,"Cup",bnd_j,"Ndn",bnd_i;
  ampo += -eta*eta*eta*a3,"Cdagup",bnd_i,"Cup",bnd_j,"Ndn",bnd_j;
  ampo += eta*eta*eta*a3/2,"Cdagup",bnd_i,"Cup",bnd_j;
  ampo += eta*eta*eta*2*a3,"Cdagup",bnd_j,"Cup",bnd_i,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += -eta*eta*eta*a3,"Cdagup",bnd_j,"Cup",bnd_i,"Ndn",bnd_i;
  ampo += -eta*eta*eta*a3,"Cdagup",bnd_j,"Cup",bnd_i,"Ndn",bnd_j;
  ampo += eta*eta*eta*a3/2,"Cdagup",bnd_j,"Cup",bnd_i;
  ampo += eta*eta*eta*2*a3,"Cdagdn",bnd_i,"Cdn",bnd_j,"Nup",bnd_i,"Nup",bnd_j;
  ampo += -eta*eta*eta*a3,"Cdagdn",bnd_i,"Cdn",bnd_j,"Nup",bnd_i;
  ampo += -eta*eta*eta*a3,"Cdagdn",bnd_i,"Cdn",bnd_j,"Nup",bnd_j;
  ampo += eta*eta*eta*a3/2,"Cdagdn",bnd_i,"Cdn",bnd_j;
  ampo += eta*eta*eta*2*a3,"Cdagdn",bnd_j,"Cdn",bnd_i,"Nup",bnd_i,"Nup",bnd_j;
  ampo += -eta*eta*eta*a3,"Cdagdn",bnd_j,"Cdn",bnd_i,"Nup",bnd_i;
  ampo += -eta*eta*eta*a3,"Cdagdn",bnd_j,"Cdn",bnd_i,"Nup",bnd_j;
  ampo += eta*eta*eta*a3/2,"Cdagdn",bnd_j,"Cdn",bnd_i;
  ampo += -eta*eta*eta*eta*4*a4,"Nup",bnd_i,"Nup",bnd_j,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += eta*eta*eta*eta*2*a4,"Nup",bnd_i,"Nup",bnd_j,"Ndn",bnd_i;
  ampo += eta*eta*eta*eta*2*a4,"Nup",bnd_i,"Nup",bnd_j,"Ndn",bnd_j;
  ampo += -eta*eta*eta*eta*a4,"Nup",bnd_i,"Nup",bnd_j;
  ampo += eta*eta*eta*eta*2*a4,"Nup",bnd_i,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += -eta*eta*eta*eta*a4,"Nup",bnd_i,"Ndn",bnd_i;
  ampo += -eta*eta*eta*eta*a4,"Nup",bnd_i,"Ndn",bnd_j;
  ampo += eta*eta*eta*eta*a4/2,"Nup",bnd_i;
  ampo += eta*eta*eta*eta*2*a4,"Nup",bnd_j,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += -eta*eta*eta*eta*a4,"Nup",bnd_j,"Ndn",bnd_i;
  ampo += -eta*eta*eta*eta*a4,"Nup",bnd_j,"Ndn",bnd_j;
  ampo += eta*eta*eta*eta*a4/2,"Nup",bnd_j;
  ampo += -eta*eta*eta*eta*a4,"Ndn",bnd_i,"Ndn",bnd_j;
  ampo += eta*eta*eta*eta*a4/2,"Ndn",bnd_i;
  ampo += eta*eta*eta*eta*a4/2,"Ndn",bnd_j;
  ampo += -eta*eta*eta*eta*a4/4,"Id",bnd_i;
  
  return ampo;
}
