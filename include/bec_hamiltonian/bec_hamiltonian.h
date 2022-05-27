//
// Created by Thomas Prosser on 27.05.22.
//

#include "solver.h"

#ifndef BEC_LATTICE_KICK_BEC_HAMILTONIAN_H
#define BEC_LATTICE_KICK_BEC_HAMILTONIAN_H

namespace bec {

  template<class T, class M, class V, bool Lattice, bool Trap>
  class Hamiltonian{
  public:
    Hamiltonian(const int& N, const T& L, const T& g, const V& psi, int upper_limit=1000, bool withit=false);
    V find_groundstate();

    template<class In>
    bool check(In first, In last, In newvec);

  private:
    Solver<T, M, V, Lattice, Trap> s;
    V output_, result_;
    bool wi_;
    int ul_, it_;
  };

}

#endif // BEC_LATTICE_KICK_BEC_HAMILTONIAN_H
