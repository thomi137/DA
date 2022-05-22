#ifndef BEC_H
#define BEC_H

#include<cmath>
#include<algorithm>
#include<functional>
#include<vector>
#include<iterator>
#include<complex>
#include<iostream>
#include "fftw3.h"

//a handy definition...
const double pi = 4*atan(1.0);

	const std::complex<double> I(0,1);
	

	template<class In, class T>
    	inline void mod_square(In first, In last, T& init)
    	{	    
	typedef typename In::value_type value_type;
	while(first != last){
		init+=std::norm<value_type>(*first);
		++first;
 		}
    	}
  
   	template<class In, class Out>    
   	inline void conj_vector(In first, In last, Out res){
		typedef typename In::value_type value_type;
   	 	while(first != last){
	    		*res++ = std::conj<value_type>(*first);
  			++first;
		}
  	}



namespace bec{
	
	template<class T, bool Lattice, bool Trap>
	struct Potential{
		public:
			Potential();
			Potential(T x):x_(x), q_(0.){}
			Potential(T x, T q):x_(x), q_(q){}
			inline T operator()(){}
			
		private:
			T x_, q_;
	};	
	
	template<class T>
	struct Potential<T, true, true>{
		inline T operator()(T x, T q = 0.){
			T sinx = sin(q*x);
			return 0.5*x*x+0.5*q*q*sinx*sinx;
		}
	};
	
	template<class T>
	struct Potential<T, false, true>{
		inline T operator()(T x, T q = 0.){
			return 0.5*x*x;
		}
	};
                
	template<class T>
        struct Potential<T, true, false>{
		inline T operator()(T x, T q = 0.){
                    T sinx = sin(q*x);
                    return 0.5*q*q*sinx*sinx;
		}
            };
        
        
	template<class T>
            struct Potential<T, false, false>{
		inline T operator()(T x, T q = 0.){
                    return 0.;
		}
            };
        
             
 	template<class VECTOR, bool Lattice, bool Trap>
	double energy(VECTOR psi, Potential<double, Lattice, Trap> pot, double L, double coupling_strength = 0.){
		int N = psi.size();
		typedef typename VECTOR::value_type value_type;
		double h = L/double(N);
 		value_type integral = 0.;
		integral+=(-0.5*conj(psi[0])*(psi[1]-2.*psi[0]+psi[N-1])/(h*h)+
				0.5*coupling_strength*conj(psi[0])*norm(psi[0])*(psi[0])+pot(L/2.)*norm(psi[0]));
		integral+=(-0.5*conj(psi[N-1])*(psi[0]-2.*psi[N-1]+psi[N-2])/(h*h)+
			0.5*coupling_strength*conj(psi[N-1])*norm(psi[N-1])*(psi[N-1])+pot(L/2.-double(N-1)*L/double(N))*norm(psi[N-1]) );
		for(int i = 1; i<N-1; ++i){
			double xpos = L/2.-double(i)*L/double(N);
			integral+=(3.-((i%2==0)?1.:-1.))*(-0.5*conj(psi[i])*(psi[i+1]-2.*psi[i]+psi[i-1])/(h*h)+
				0.5*coupling_strength*conj(psi[i])*norm(psi[i])*psi[i]+pot(xpos)*norm(psi[i]));
		}		
   	  	return real(L/(3.*double(N))*integral);

        }

}//namespace bec


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
template<class LTYPE>
Momentum_Operator<LTYPE>::Momentum_Operator(LTYPE length):length_(length){}

template<class LTYPE>
 template<class VECTOR>
  inline VECTOR Momentum_Operator<LTYPE>::operator*(VECTOR psi) const{
	for(int i = 0; i<psi.size(); ++i)
		psi[i] *= (double(i)<double(psi.size())/2.)?2.*pi/(length_)*double(i): -2.*pi/(length_)*double(psi.size()-i);
	return psi;
}
namespace Mean_Op{



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
		value_type secder = -1./(h*h)*(psi[(i<psi.size()-1)?i+1:0]-2.*psi[i]+psi[(i==0)?psi.size()-1:i-1]);
		integral+=(3.-((i%2==0)?1.:-1.))*conj(psi[i])*secder;
	}
	integral*=length/(3.*double(psi.size() ) );
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
	inline typename VECTOR::value_type operator()(VECTOR psi) const;
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


}//namespace Mean_Op
#endif    
    
    
  
