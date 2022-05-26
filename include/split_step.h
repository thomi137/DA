#ifndef SPLIT_STEP_H
#define SPLIT_STEP_H

#include <complex>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>

#include "fft.h"
#include "bec.h"

namespace split_step{

	template<class T, bool Lattice, bool Trap>
	inline std::vector<std::complex<double>> S3 (const T L,
                                               const T deltat,
                                               const T g, std::vector<std::complex<double>> data,
                                               bec::Potential<T, Lattice, Trap> pot,
                                               const T omega=1.) {
    double q = bec::pi*40./L;
    std::complex<double>I(0,1);
    unsigned int N = data.size();
    std::vector<std::complex<double> > kvec(N);
    std::vector<std::complex<double> > kfac(N);
                      std::vector<std::complex<double> > transform(N);
    for (unsigned i = 0; i<N; ++i) {
      kvec[i]= (double(i)<double(N)/2.)?2.*pi/(L)*double(i): -2.*pi/(L)*double(N-i);
              kfac[i]= exp(-bec::I*deltat*omega/2.*kvec[i]*kvec[i]/2.);
    }
    transform=FourierTransform<fft::forward, fft::Estimate>(data);
    std::transform(transform.begin(), transform.end(), kfac.begin(), transform.begin(), std::multiplies<std::complex<double> >() );
    std::transform(transform.begin(), transform.end(), transform.begin(),std::bind2nd(std::divides<std::complex<double> >(), N) );
    data=FourierTransform<fft::backward, fft::Estimate>(transform);

    for(unsigned int i = 0; i<N; ++i){
        double xpos = L*0.5-double(i)*L/double(N);
        double sinx = std::sin(q*xpos);
        data[i]*=exp(-I*omega*deltat*(g*std::norm(data[i])+pot(xpos, q) ) );
    }

    transform=FourierTransform<fft::forward, fft::Estimate>(data);
    std::transform(transform.begin(), transform.end(), kfac.begin(), transform.begin(), std::multiplies<std::complex<double> >() );
    std::transform(transform.begin(), transform.end(), transform.begin(),bind2nd(std::divides<std::complex<double> >(), N) );
    data=FourierTransform<fft::backward, fft::Estimate>(transform);

    return data;
		}

    template<class T, bool Lattice, bool Trap>
    inline std::vector<std::complex<double> > S3imag (const T L,
                                                      const T deltat,
                                                      const T g, std::vector<std::complex<double>> data,
                                                      bec::Potential<T, Lattice, Trap> pot) {
      const double pi = 4.*atan(1.);
      double q = pi*40./L;
      std::complex<double>I(0,1);
      unsigned int N = data.size();
      std::vector<std::complex<double> > kvec(N);
      std::vector<std::complex<double> > kfac(N);
      std::vector<std::complex<double> > transform(N);
      for (unsigned i = 0; i<N; ++i) {
          kvec[i]= (double(i)<double(N)/2.)?2.*pi/(L)*double(i): -2.*pi/(L)*double(N-i);
          kfac[i]= exp(-deltat/2.*kvec[i]*kvec[i]/2.);
      }
      transform=FourierTransform<fft::forward, fft::Estimate>(data);
      std::transform(transform.begin(), transform.end(), kfac.begin(), transform.begin(), std::multiplies<std::complex<double> >() );
      std::transform(transform.begin(), transform.end(), transform.begin(),std::bind2nd(std::divides<std::complex<double> >(), N) );
      data=FourierTransform<fft::backward, fft::Estimate>(transform);

      for(unsigned int i = 0; i<N; ++i){
          double xpos =L*0.5-double(i)*L/double(N);
          double sinx =sin(q*xpos);
          data[i]*=exp(-deltat*(g*std::norm(data[i])+pot(xpos, q) ) );
      }

      transform=FourierTransform<fft::forward, fft::Estimate>(data);
      std::transform(transform.begin(), transform.end(), kfac.begin(), transform.begin(), std::multiplies<std::complex<double> >() );
      std::transform(transform.begin(), transform.end(), transform.begin(),bind2nd(std::divides<std::complex<double> >(), N) );
      data=FourierTransform<fft::backward, fft::Estimate>(transform);

      return data;
  }

	template<class T, bool Lattice, bool Trap>	
	inline std::vector<std::complex<double> > S7 (const T L, const T deltat, const T g, std::vector<std::complex<double>
	> data, bec::Potential<T, Lattice, Trap> pot) {
    //Bandrauk and Shen, J. Phys. A, 27:7147
    const double omega[4] =
            {1.3151861047504085, -1.1776799841887, 0.2355733213359357, 0.784513610477560};
    data=S3(L, deltat, g, data, pot, omega[3]);
    data=S3(L, deltat, g, data, pot, omega[2]);
    data=S3(L, deltat, g, data, pot, omega[1]);
    data=S3(L, deltat, g, data, pot, omega[0]);
    data=S3(L, deltat, g, data, pot, omega[1]);
    data=S3(L, deltat, g, data, pot, omega[2]);
    data=S3(L, deltat, g, data, pot, omega[3]);
    return data;
  }
} //namespace split_step

#endif	
