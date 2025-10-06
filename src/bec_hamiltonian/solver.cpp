//
// Created by Thomas Prosser on 27.05.22.
//

#include "solver.h"


namespace bec {
  /**
   * Compute eigenvalues of matrices. For documentation, have a look at:
   * http://www.netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_gaaa6df51cfd92c4ab08d41a54bf05c3ab.html
   */
  extern "C" void dstev_(char* jobz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info);

  Solver::Solver(const int &N, const double &L, const double &g, const Vector &psi)
      :N_(N), L_(L), g_(g), ldz_(N_), jobz_('V'), deltax(L_/double(N_) ), q_(40.*pi/L_), psi_(psi), d_(N_), e_(N_), Rmat_(N, N){
    work_ = new double[2*N_-2];
    init(bec::Potential<bec_t , true, true>());
  }

  void Solver::init(bec::Potential<bec_t ,true, true> pot){
    for(size_t i=0; i<Rmat_.size1();++i){
      double xpos =L_*0.5-double(i)*L_/double(N_);
      d_[i]=1/(deltax*deltax) + g_*psi_[i]*psi_[i] + pot(xpos, q_);
      e_[i]=-0.5/(deltax*deltax);
    }
    e_[N_-1]=0;
  }

  void Solver::diag(Vector &res){
    //the work gets done here!!!
    dstev_(&jobz_, &N_, &d_[0], &e_[0], &Rmat_(0,0), &ldz_, work_, &info_);
    //eigenvectors are now stored in Rmat_, we need the lowest one and extract it...
    double* point =  &Rmat_(0,0);
    copy(point, point+N_ , res.begin());
  }

  void Solver::re_init(const Vector & start, bec::Potential<bec_t , true, true>& pot){
    psi_.clear();
    copy(start.begin(), start.end(), back_inserter(psi_));
    init(pot);
  }

}