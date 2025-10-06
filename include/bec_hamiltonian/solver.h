//
// Created by Thomas Prosser on 27.05.22.
//

#ifndef BEC_LATTICE_KICK_SOLVER_H
#define BEC_LATTICE_KICK_SOLVER_H

#include "bec.h"

namespace bec {

  class Solver {
  public:
    Solver(const int& N, const double& L, const double& g, const Vector &psi);
    ~Solver() = default;

    /**
       *
       * @tparam Denotes the type of the output
       * @param res output vector
     */
    void diag(Vector &res);

    void init(bec::Potential<bec_t , true, true> pot);
    void re_init(const Vector & start, bec::Potential<bec_t , true, true> &pot);

  private:
    int N_, ldz_, lwork_, info_;
    bec_t L_, g_, deltax, q_;
    Vector psi_;
    double* work_;
    Vector d_, e_;
    char jobz_;
    Matrix Rmat_;
  };

}; // namespace bec

#endif // BEC_LATTICE_KICK_SOLVER_H
