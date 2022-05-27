//
// Created by Thomas Prosser on 27.05.22.
//

#include "solver.h"

/**
 * Compute eigenvalues of matrices. For documentation, have a look at:
 * http://www.netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_gaaa6df51cfd92c4ab08d41a54bf05c3ab.html
 */
extern "C" void dstev_(char* jobz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info);

namespace bec {

  template<class T, class M, class V, bool Lattice, bool Trap>
  Solver<T, M, V, Lattice, Trap>::Solver(const int& N, const double&L, const double& g, const V& psi)
      :N_(N), L_(L), g_(g), ldz_(N_), jobz_('V'), deltax(L_/double(N_) ), q_(40.*pi/L_), psi_(psi), d_(N_), e_(N_), Rmat_(N, N){
    work_ = std::unique_ptr<T[]>(new T[2*N_-2]);
    init(bec::Potential<T, Lattice, Trap>());
  }

  template<class T, class M, class V, bool Lattice, bool Trap>
  void Solver<T, M, V, Lattice, Trap>::init(bec::Potential<T, Lattice, Trap> pot){
    for(size_t i=0; i<Rmat_.size1();++i){
      double xpos =L_*0.5-double(i)*L_/double(N_);
      d_[i]=1/(deltax*deltax) + g_*psi_[i]*psi_[i] + pot(xpos, q_);
      e_[i]=-0.5/(deltax*deltax);
    }
    e_[N_-1]=0;
  }

  template<class T, class M, class V, bool Lattice, bool Trap>
  template<class Out>
  void Solver<T, M, V, Lattice, Trap>::diag(Out res){
    //the work gets done here!!!
    dstev_(&jobz_, &N_, &d_[0], &e_[0], &Rmat_(0,0), &ldz_, work_, &info_);
    //eigenvectors are now stored in Rmat_, we need the lowest one and extract it...
    double* point = std::make_unique<double>( &Rmat_(0,0) );
    copy(point, point+N_, res);
  }

}