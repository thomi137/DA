//
// Created by Thomas Prosser on 24.05.22.
//

#ifndef MEAN_OP_H
#define MEAN_OP_H

#include<cmath>
#include<algorithm>
#include<functional>
#include<vector>
#include<iterator>
#include<complex>
#include<iostream>

#include "fftw3.h"
#include "./bec.h"

namespace bec {

  template<class VECTOR>
  inline typename VECTOR::value_type mean_momentum(VECTOR psi, double length = 10.){
    typedef typename VECTOR::value_type value_type;
    value_type integral = 0;
    double h = length/double(psi.size() );
    for (int i = 0; i<psi.size();++i){
      value_type derivative = -1.*I/(2.*h)*(psi[(i<psi.size()-1)?i+1:0]-psi[(i==0)?psi.size()-1:i-1]);
      integral+=(3.-((i%2==0)?1.:-1.))*std::conj(psi[i])*derivative;
    }
    integral*=length/(3.*double(psi.size() ) );
    return integral;
  }

  template<class VECTOR>
  inline typename VECTOR::value_type mean_square_momentum(VECTOR psi, double length = 10.){
    typedef typename VECTOR::value_type value_type;
    value_type integral = 0;
    double h = length/double(psi.size() );
    for (int i = 0; i<psi.size();++i){
      value_type secder = -1./(h*h) * (psi[(i<psi.size()-1)?i+1:0]-2. * psi[i]+psi[(i==0)?psi.size()-1:i-1]);
      integral+=(3.-((i%2==0)?1.:-1.)) * std::conj(psi[i])*secder;
    }
    integral*=length/(3.*double(psi.size()));
    return integral;
  }

  template<class VECTOR>
  inline typename VECTOR::value_type delta_momentum(VECTOR psi, double length = 10.){
    typedef typename VECTOR::value_type value_type;
    value_type xsq = mean_momentum(psi, length);
    return (mean_square_momentum(psi, length)-xsq*xsq);
  }

  template<class VECTOR>
  inline double mean_position(VECTOR psi, double length = 10.){
    double integral = 0.;
    for(int i = 0; i < psi.size(); ++i){
      double xpos = length/2.-double(i)* length/double(psi.size());
      integral+=(3.-((i%2==0)?1.:-1.))*xpos*std::norm(psi[i]);
    }
    integral*=length/(3.*double(psi.size() ) );
    return integral;
  }

  struct Momentum{};
  struct Position{};

  template<class OP, class VECTOR>
  struct Operator_Mean{
    inline typename VECTOR::value_type operator() (VECTOR psi) const;
  };

  template<class VECTOR>
  inline double mean_square_position(VECTOR psi, double length = 10.){
    double integral = 0.;
    for(int i = 0; i < psi.size(); ++i){
      double xpos = length/2.-double(i)* length/double(psi.size());
      integral+=(3.-((i%2==0)?1.:-1.))*xpos*xpos*std::norm(psi[i]);
    }
    integral*=length/(3.*double(psi.size() ) );
    return integral;
  }

  template<class VECTOR>
  inline double delta_position(VECTOR psi, double length = 10.){
    typedef typename VECTOR::value_type value_type;
    value_type xsq = mean_position(psi, length);
    return (mean_square_position(psi, length)-std::abs(xsq*xsq));
  }

} //namespace bec

#endif // MEAN_OP_H
