#include "itensor/all.h"
#include "SpinChain.h"

using namespace itensor;

MPO Hamiltonian::SpinChain()
{
  auto ampo = AutoMPO(sites);
  for(int b = 1; b < N + periodic; ++b)
    {
      ampo += 0.5*J1,"S+",b,"S-",b%N+1;
      ampo += 0.5*J1,"S-",b,"S+",b%N+1;
      ampo += J1,"Sz",b,"Sz",b%N+1;
    }
  for(int b = 1; b < N-1 + 2*periodic; ++b)
    {
      ampo += 0.5*J2,"S+",b,"S-",(b+1)%N+1;
      ampo += 0.5*J2,"S-",b,"S+",(b+1)%N+1;
      ampo += J2,"Sz",b,"Sz",(b+1)%N+1;
    }
  auto H = toMPO(ampo,{"Cutoff=",cutoff});

  return H;
}
