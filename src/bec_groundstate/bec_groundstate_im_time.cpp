#include<iostream>
#include<algorithm>
#include<iterator>
#include<cmath>
#include<complex>

#include "bec_groundstate.h"
#include"fftw3.h"
#include"bec.h"
#include"split_step.h"

using namespace boost::numeric::ublas;

namespace bec {
  void BecGroundstateImTime::execute() {

    bec_t mynorm = 0.;

    BecParameters params;
    auto [period, N, K, g, maxt, L] = params;
    bec_t castN = boost::numeric_cast<bec_t>(N);
    bec_t  castPeriod = boost::numeric_cast<bec_t>(period);



    const double deltat = 2. * pi / castPeriod;

    Cvector groundstate(N);
    for (int i = 0; i < groundstate.size(); ++i) {
      double xpos = L * 0.5 - double(i) * L / castN;
      groundstate[i] = 1. / sqrt(sqrt(pi)) * exp(-xpos * xpos / 2.);
    }

    for (int i = 0; i < 20000; ++i) {

      // Calculate S3
      groundstate = split_step::S3imag(L, deltat, g, groundstate,
                                       bec::Potential<double, true, false>());

      // Normalize groundstate.
      // std::bind2nd is deprecated since C++11 so will go for lambdas. Looks better anyway.
      mynorm = L_2_norm(groundstate.begin(), groundstate.end(), L);
      transform(groundstate.begin(), groundstate.end(), groundstate.begin(),
                [&mynorm](complex_type el) -> complex_type { return el / mynorm; });
    }

    Cvector data(N);
    Cvector mom(N);
    copy(groundstate.begin(), groundstate.end(), data.begin());

    std::cout << "System size" << '\t' << "Points" << '\t' << "Coupling" << '\t'
              << "deltat" << '\t' << "Kick strength" << '\t' << "Max time" << '\t'
              << "Period" << '\n';
    std::cout << L << '\t' << N << '\t' << g << '\t' << deltat << '\t' << K
              << '\t' << maxt << '\t' << period << '\t';

    // Kicking!
    Cvector kick(N);
    for (size_t i = 0; i < N; ++i) {
      bec_t xpos = L * 0.5 - double(i) * L / castN;
      kick[i] = exp(-I * K * sin(xpos));
    }

    transform(data.begin(), data.end(), kick.begin(), data.begin(),
              std::multiplies<complex_type>());

    double t = 0;
    int stepcounter = 0;

    for (; t < maxt; t += deltat, ++stepcounter) {

      data = split_step::S3(L, deltat, g, data, Potential<double, true, false>());
      if (!(stepcounter % (period))) {

        mom = FourierTransform<fft::forward, fft::Estimate>(data);
        for (int i = 0; i < N; ++i) {
          int index = i < N / 2 ? i + N / 2. : i - N / 2.;
          std::cout << norm(data[i]) << "\t" << norm(mom[index]) << '\n';
        }
      }
    }
  }
} // namespace bec
