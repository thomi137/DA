//
// Created by Thomas Prosser on 27.05.22.
//

#ifndef BEC_LATTICE_KICK_SOLVER_H
#define BEC_LATTICE_KICK_SOLVER_H

#include "bec.h"

namespace bec {

template<class T, class M, class V, bool Lattice, bool Trap>

  class Solver {
  public:
    Solver(const int& N, const double& L, const double& g, const V& psi);
    ~Solver() = default;

    /**
       *
       * @tparam Denotes the type of the output
       * @param res output vector
     */
    template<class Out>
    void diag(Out res);
    void init(bec::Potential<T, Lattice, Trap> pot);

    void re_init(const V& start, bec::Potential<T, Lattice, Trap> pot);

  private:
    int N_, ldz_, lwork_, info_;
    T L_, g_, deltax, q_;
    V psi_;
    T* work_;
    V d_, e_;
    char jobz_;
    V Rmat_;
  };

}; // namespace bec

#endif // BEC_LATTICE_KICK_SOLVER_H
