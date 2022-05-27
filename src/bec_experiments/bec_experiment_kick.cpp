//
// Created by Thomas Prosser on 27.05.22.
//

#include "bec.h"
#include "fft.h"
#include "bec_hamiltonian.h"
#include "bec_groundstate.h"

namespace bec {

  void BecExperimentKick::execute() {

    BecParameters params;
    auto [period, N, K, g, maxt, L] = params;
    bec_t castN = boost::numeric_cast<bec_t>(N);
    bec_t  castPeriod = boost::numeric_cast<bec_t>(period);

    const double deltat = 2.*pi/double(period);

    Vector groundstate(N, 100.);
    Cvector data(N);
    Hamiltonian<bec_t, Matrix, Vector, true, true> GP(N, L, g, groundstate, 50000, true);
    groundstate=GP.find_groundstate();
    copy(groundstate.begin(), groundstate.end(), data.begin());
    std::cout<<"System size"<<'\t'<<"Points"<<'\t'<<"Coupling"<<'\t'<<"deltat"<<std::endl;
    std::cout<<L<<'\t'<<N<<'\t'<<g<<'\t'<<deltat<<std::endl;

    Cvector mom(N);
    mom= FourierTransform<fft::forward,fft::Estimate>(data);
    for(int i = 0 ; i < mom.size(); ++i)
      std::cout<<norm(mom[i])<<'\n';
  }


}