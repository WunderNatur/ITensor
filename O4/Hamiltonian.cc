#include "O4.h"

using namespace itensor;

MPO O4::O4Hamiltonian()
{
  auto ampoO4 = AutoMPO(sites);
  for(int i = 1; i <= N; ++i)
    {
      ampoO4 += U,"Nupdn",i;
      ampoO4 += -U/2,"Nup",i;
      ampoO4 += -U/2,"Ndn",i;
      ampoO4 += U/4,"Id",i;
    }
  for(int b = 1; b < N + periodic; ++b)
    {
      ampoO4 += -0.25,"Id",b;
      ampoO4 += -a1/2,"Cdagup",b,"Cup",b%N+1;
      ampoO4 += -a1/2,"Cdagup",b%N+1,"Cup",b;
      ampoO4 += -a1/2,"Cdagdn",b,"Cdn",b%N+1;
      ampoO4 += -a1/2,"Cdagdn",b%N+1,"Cdn",b;
      ampoO4 += -a2,"Cdagup",b,"Cup",b%N+1,"Cdagdn",b,"Cdn",b%N+1;
      ampoO4 += -a2,"Cdagup",b,"Cup",b%N+1,"Cdagdn",b%N+1,"Cdn",b;
      ampoO4 += -a2,"Cdagup",b%N+1,"Cup",b,"Cdagdn",b,"Cdn",b%N+1;
      ampoO4 += -a2,"Cdagup",b%N+1,"Cup",b,"Cdagdn",b%N+1,"Cdn",b;
      ampoO4 += a2,"Nup",b,"Nup",b%N+1;
      ampoO4 += -a2/2,"Nup",b;
      ampoO4 += -a2/2,"Nup",b%N+1;
      ampoO4 += a2/4,"Id",b;
      ampoO4 += a2,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += -a2/2,"Ndn",b;
      ampoO4 += -a2/2,"Ndn",b%N+1;
      ampoO4 += a2/4,"Id",b;
      ampoO4 += 2*a3,"Cdagup",b,"Cup",b%N+1,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += -a3,"Cdagup",b,"Cup",b%N+1,"Ndn",b;
      ampoO4 += -a3,"Cdagup",b,"Cup",b%N+1,"Ndn",b%N+1;
      ampoO4 += a3/2,"Cdagup",b,"Cup",b%N+1;
      ampoO4 += 2*a3,"Cdagup",b%N+1,"Cup",b,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += -a3,"Cdagup",b%N+1,"Cup",b,"Ndn",b;
      ampoO4 += -a3,"Cdagup",b%N+1,"Cup",b,"Ndn",b%N+1;
      ampoO4 += a3/2,"Cdagup",b%N+1,"Cup",b;
      ampoO4 += 2*a3,"Cdagdn",b,"Cdn",b%N+1,"Nup",b,"Nup",b%N+1;
      ampoO4 += -a3,"Cdagdn",b,"Cdn",b%N+1,"Nup",b;
      ampoO4 += -a3,"Cdagdn",b,"Cdn",b%N+1,"Nup",b%N+1;
      ampoO4 += a3/2,"Cdagdn",b,"Cdn",b%N+1;
      ampoO4 += 2*a3,"Cdagdn",b%N+1,"Cdn",b,"Nup",b,"Nup",b%N+1;
      ampoO4 += -a3,"Cdagdn",b%N+1,"Cdn",b,"Nup",b;
      ampoO4 += -a3,"Cdagdn",b%N+1,"Cdn",b,"Nup",b%N+1;
      ampoO4 += a3/2,"Cdagdn",b%N+1,"Cdn",b;
      ampoO4 += -4*a4,"Nup",b,"Nup",b%N+1,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += 2*a4,"Nup",b,"Nup",b%N+1,"Ndn",b;
      ampoO4 += 2*a4,"Nup",b,"Nup",b%N+1,"Ndn",b%N+1;
      ampoO4 += -a4,"Nup",b,"Nup",b%N+1;
      ampoO4 += 2*a4,"Nup",b,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += -a4,"Nup",b,"Ndn",b;
      ampoO4 += -a4,"Nup",b,"Ndn",b%N+1;
      ampoO4 += a4/2,"Nup",b;
      ampoO4 += 2*a4,"Nup",b%N+1,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += -a4,"Nup",b%N+1,"Ndn",b;
      ampoO4 += -a4,"Nup",b%N+1,"Ndn",b%N+1;
      ampoO4 += a4/2,"Nup",b%N+1;
      ampoO4 += -a4,"Ndn",b,"Ndn",b%N+1;
      ampoO4 += a4/2,"Ndn",b;
      ampoO4 += a4/2,"Ndn",b%N+1;
      ampoO4 += -a4/4,"Id",b;
    }
  auto H = toMPO(ampoO4,{"Cutoff=",cutoff});

  return H;
}
