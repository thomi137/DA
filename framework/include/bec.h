#ifndef BEC_H
#define BEC_H

#include<cmath>
#include<algorithm>
#include<functional>
#include<vector>
#include<iterator>
#include<complex>
#include<iostream>
#include <boost/math/constants/constants.hpp>

#include "fftw3.h"


namespace bec {

  // We're living in a 64 bit world now. It's no longer 2003, long double is here...
  typedef double bec_t;

  // Goes without words
  const bec_t pi = boost::math::constants::pi<bec_t>();

  // Imaginary unit
  const std::complex<bec_t> I(0,1);

  /**
   * Calculates sum_i(x_i*x_i^star), or the squared modulus of a vector.
   *
   * @tparam In
   * @tparam T  Usually a double or complex, but could change.
   * @param first Iterator on the input vectors start
   * @param last Iterator on the input vectors end
   * @param init Allocated result
   */
  template<class In, class T>
  inline void mod_square(In first, In last, T& init) {
    typedef typename In::value_type value_type;
    while(first != last){
      init+=std::norm<value_type>(*first);
      ++first;
    }
  }

  /**
   * Requires more documentation. In place operations
   *
   * @tparam In
   * @tparam Out Reusult vector. TODO: Make consistent with mod_square
   * @tparam T Usually double or complex, but could change
   * @param first Iterator on the input vectors start
   * @param last Iterator on the input vectors end
   */
  template<class In, class Out>
  inline void conj_vector(In first, In last, Out res) {
    typedef typename In::value_type value_type;
    while(first != last){
      *res++ = std::conj<value_type>(*first);
      ++first;
    }
  }

  /**
   * Template class definition for the BEC potential as used in the thesis (cf. chapter 1)
   *
   * @tparam T numeric type.
   * @tparam Lattice harmonic lattice (sine-shaped)
   * @tparam Trap Whether the system is trapped in a harmonic oscillator potential.
   *
   * q denotes the wave vector of the incident laser in an optical trap. Cf. eq. 1.5 in the thesis for more
   * information.
   */
  template<class T, bool Lattice, bool Trap>
  struct Potential{
    public:
      Potential();
      Potential(T x):x_(x), q_(0.){}
      Potential(T x, T q):x_(x), q_(q){}
      inline T operator() () {}

    private:
      T x_, q_;
  };

	template<class T>
	struct Potential<T, true, true>{
		inline T operator()(T x, T q){
			T sinx = sin(q*x);
			return (0.5*x*x) + (0.5*q*q*sinx*sinx);
		}
	};

  // Trapped potential without a lattice
	template<class T>
	struct Potential<T, false, true>{
		inline T operator()(T x, T q = 0.){
			return 0.5*x*x;
		}
	};

  // Free potential in a lattice
	template<class T>
	struct Potential<T, true, false> {
    inline T operator()(T x, T q = 0.){
      T sinx = sin(q*x);
      return 0.5*q*q*sinx*sinx;
		}
  };
        
  // Free particle. No potential.
	template<class T>
  struct Potential<T, false, false>{
    inline T operator()(T x, T q = 0.){
            return 0.;
    }
  };

 	template<class VECTOR, bool Lattice, bool Trap>
	double energy(VECTOR psi, Potential<double, Lattice, Trap> pot, double L, double coupling_strength = 0.) {
		typedef typename VECTOR::value_type value_type;
    value_type integral = 0.;

    int N = psi.size();
		double h = L/double(N);

    integral += (-0.5*conj(psi[0])*(psi[1]-2.*psi[0]+psi[N-1])/(h*h)+
				0.5*coupling_strength*std::conj(psi[0])*std::norm(psi[0])*(psi[0])+pot(L/2.)*std::norm(psi[0]));
		integral+=(-0.5*std::conj(psi[N-1]) * (psi[0]-2.*psi[N-1]+psi[N-2])/(h*h)+
			0.5*coupling_strength * std::conj(psi[N-1])* std::norm(psi[N-1])*(psi[N-1])+pot(L/2.-double(N-1)*L/double(N))* std::norm(psi[N-1]) );
		for(int i = 1; i<N-1; ++i){
			double xpos = L/2.-double(i)*L/double(N);
			integral+=(3.-((i%2==0)?1.:-1.))*(-0.5*conj(psi[i])*(psi[i+1]-2.*psi[i]+psi[i-1])/(h*h)+
				0.5*coupling_strength*conj(psi[i])*norm(psi[i])*psi[i]+pot(xpos)*norm(psi[i]));
		  }
    return real(L/(3.*double(N))*integral);
  }



  // Again... I was very young... Had to show off pointer arithmetics. Dangerous stuff not made for mortal men ;-)
  template <class In, class T>
  inline T  L_2_norm( In first, In last, T range){
    T init = 0.;
            int N = last - first;
    while(first != last)
      init+=norm(*first++);
    return sqrt(range/double(N)*init);
  }

  template<class LTYPE>
  class Momentum_Operator{
    public:
      Momentum_Operator(LTYPE length);
      template<class VECTOR>
      inline VECTOR operator*(VECTOR psi) const;

    private:
      LTYPE length_;

  };

  // Implementation of Momentum_Operator
  template<class LTYPE>
  Momentum_Operator<LTYPE>::Momentum_Operator(LTYPE length):length_(length){}

  template<class LTYPE>
   template<class VECTOR>
    inline VECTOR Momentum_Operator<LTYPE>::operator*(VECTOR psi) const{
    for(int i = 0; i<psi.size(); ++i)
      psi[i] *= (double(i) < double(psi.size())/2.)
              ? 2.*pi/(length_)*double(i)
              : -2.*pi/(length_)*double(psi.size()-i);
    return psi;
  }
} //namespace bec


#endif // BEC_H
