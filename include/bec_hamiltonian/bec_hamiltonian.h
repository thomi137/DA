//
// Created by Thomas Prosser on 27.05.22.
//

#include "solver.h"

#ifndef BEC_LATTICE_KICK_BEC_HAMILTONIAN_H
#define BEC_LATTICE_KICK_BEC_HAMILTONIAN_H

namespace bec {

  class Hamiltonian {
  public:
    Hamiltonian(const int &N, const bec_t &L, const bec_t &g, const Vector &psi, bool with_iterations);
    Vector find_groundstate();
  private:
    template<class In>
    bool check(In first, In last, In newvec);

    Solver s;
    bec::Potential<bec_t , true, true> pot;
    Vector output_, result_;
    bool wi_ = false;
    int ul_, it_ = 0;
  };

}

#endif // BEC_LATTICE_KICK_BEC_HAMILTONIAN_H
