//
// Created by Thomas Prosser on 27.05.22.
//

#include "bec_hamiltonian.h"

namespace bec {

  Hamiltonian::Hamiltonian(const int &N, const bec_t &L, const bec_t &g, const Vector &psi, bool with_iterations)
  :output_(psi), result_(N), s(N, L, g, psi), wi_(with_iterations), ul_(10000) {}

  // find_groundstate
  Vector Hamiltonian::find_groundstate(){
    for(int i = 0; i<ul_; ++i){
      s.diag(result_);
      transform(output_.begin(), output_.end(), result_.begin(), output_.begin(), std::minus<double>());
      // copy(result_.begin(), result_.end(), std::ostream_iterator<double>(std::cout,"\n"));
      if (check(output_.begin(), output_.end(), result_.begin())){
        std::cerr<<it_<<" iteration(s) needed"<<std::endl;
        return result_;
      }
      else {if(wi_){++it_;}
        copy(result_.begin(), result_.end(), output_.begin() );
        s.re_init(result_, pot);
      }
    }
    //if we arrive here, there is nothing more to do, we can only return the result.
    return result_;
  }


  template<class In>
  bool Hamiltonian::check(In first, In last, In newvec) {
    while (first != last) {
      if (abs(*first) <= 100. * std::numeric_limits<double>::epsilon()) {
        ++first, ++newvec;
      } else {
        return false;
      }
    }
    // if we arrive here, everything is fine
    return true;
  }

}
