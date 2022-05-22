#ifndef BEC_GROUNDSTATE_H
#define BEC_GROUNDSTATE_H

#include<iostream>
#include<cmath>
#include<vector>
#include<iterator>
#include<algorithm>
#include<limits>
#include<functional>
#include"bec.h"

//#include<ietl/lanczos.h>
//#include<ietl/interface/ublas.h>
//#include<ietl/vectorspace.h>
//#include<boost/random/lagged_fibonacci.hpp>

	 


//MATRIX should be some Matrix type provided by MTL or ublas.
using namespace std;



#ifndef IETL_LANCZOS_H
//two VERY handy declarations...
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);	
extern "C" void dstev_(char* jobz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info); 
#endif


//declarations that are of no great use as of now...
struct FullDiag {};
struct Lanczos {};


template<class T, class MATRIX, class VECTOR, class DiagPolicy, bool Lattice, bool Trap>
class solver{};

//we specialize...
template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
class solver<T, MATRIX, VECTOR, FullDiag, Lattice,Trap>{
public:
    solver(const int& N, const double& L, const double& g, const VECTOR& psi);
    ~solver();
    
    template<class Out>
    void diag(Out res);
    
    void init(bec::Potential<T, Lattice, Trap> pot);
    
    void re_init(const VECTOR& start, bec::Potential<T, Lattice, Trap> pot);
    
private:
    int N_, ldz_, lwork_, info_;
    T L_, g_, deltax, q_;
    VECTOR psi_;
    T* work_;
    VECTOR d_, e_;
    char jobz_;
    MATRIX Rmat_;
};//specialisation for the full diagonalisation
/*
//we specialize for Lanczos...
template<class MATRIX, class VECTOR>
class solver<MATRIX, VECTOR, Lanczos>{
	public:
	#include<ietl/lanczos.h>
	#include<ietl/interface/ublas.h>
	#include<ietl/vectorspace.h>
	#include<boost/random/lagged_fibonacci.hpp>

	 	solver(const int& N, const double& L, const double& g, const VECTOR& psi);
		//~solver();
		//default copy constructor seems ok (at least for now...)
		void init();
		template<class Out>
		void diag(Out res);
		void re_init(const VECTOR& start);

	private:
		typedef VECTOR::value_type value_type;
		typedef boost::numeric::ublas::vector<value_type> uvector;
		typedef boost::lagged_fibonacci607 Gen;
		typedef ietl::vectorspace<uvector> Vecspace;
		Gen mygen;
		Vecspace vec;
		ietl::Info<value_type> info;
	  	int N_;
		double L_, g_, deltax, q_;
		VECTOR psi_, eigen_, err_;
		VECTOR::iterator start, end;
		std::vector<int> mult_;
		std::vector<uvector> eigenvec_;
		MATRIX Rmat_;

};//specialisation for Lanczos.
		*/

//The little class that could...
template<class T, class MATRIX, class VECTOR, class DiagPolicy, bool Lattice, bool Trap>
class Hamiltonian{
public:
        Hamiltonian(const int& N, const T& L, const T& g, const VECTOR& psi, int upper_limit=1000, bool withit=false);
        //copy constructor not yet needed, default can be used. See Sutter on exception-safety.
        //~Hamiltonian(); 
        //default constructor ok, the only critical part is in solver, which is destroyed when Hamiltonian is...
        VECTOR find_groundstate();
        template<class In>
        bool check(In first, In last, In newvec);
    
private:
        solver<T, MATRIX, VECTOR, DiagPolicy, Lattice, Trap> s;
        VECTOR output_, result_;
        bool wi_;
        int ul_, it_;
};//base class


//////////////////////////////////////////////Implementation starts here/////////////////////////////////////////////////////////////

//solver 1:
template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
solver<T, MATRIX, VECTOR, FullDiag, Lattice, Trap>::solver(const int& N, const double&L, const double& g, const VECTOR& psi)
:N_(N), L_(L), g_(g), ldz_(N_), jobz_('V'), deltax(L_/double(N_) ), q_(40.*pi/L_), psi_(psi), d_(N_), e_(N_), Rmat_(N, N){
    work_ = new double[2*N_-2];
    init(bec::Potential<T, Lattice, Trap>());
}


template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
solver<T, MATRIX, VECTOR, FullDiag, Lattice, Trap>::~solver(){
    delete[] work_;
}

template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
void solver<T, MATRIX, VECTOR, FullDiag, Lattice, Trap>::init(bec::Potential<T, Lattice, Trap> pot){
    for(size_t i=0; i<Rmat_.size1();++i){
        double xpos =L_*0.5-double(i)*L_/double(N_);
        d_[i]=1/(deltax*deltax)+g_*psi_[i]*psi_[i]+pot(xpos, q_);
        e_[i]=-0.5/(deltax*deltax);
    }
    e_[N_-1]=0;
}

template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
 template<class Out>
  void solver<T, MATRIX, VECTOR, FullDiag, Lattice, Trap>::diag(Out res){
      //the work gets done here!!!
      dstev_(&jobz_, &N_, &d_[0], &e_[0], &Rmat_(0,0), &ldz_, work_, &info_);
      //eigenvectors are now stored in Rmat_, we need the lowest one and extract it...
      double* point = &Rmat_(0,0);
      copy(point, point+N_, res);
  }
      
template<class T, class MATRIX, class VECTOR, bool Lattice, bool Trap>
void solver<T, MATRIX, VECTOR, FullDiag, Lattice, Trap>::re_init(const VECTOR& start, bec::Potential<T, Lattice, Trap> pot){
    psi_.clear();
    copy(start.begin(), start.end(), back_inserter(psi_));
    init(pot);
}
//this concludes the implementation of solver 1.

//Hamiltonian -- constructor
template<class T, class MATRIX, class VECTOR, class DiagPolicy, bool Lattice, bool Trap>
Hamiltonian<T, MATRIX, VECTOR, DiagPolicy, Lattice, Trap>::Hamiltonian(const int& N, const T& L, const T& g, const VECTOR& psi, int upper_limit, bool
		withit)
		:wi_(withit), output_(psi), result_(N), s(N, L, g, psi), ul_(upper_limit), it_(0){}

// find_groundstate
template<class T, class MATRIX, class VECTOR, class DiagPolicy, bool Lattice, bool Trap>
VECTOR Hamiltonian<T, MATRIX, VECTOR, DiagPolicy, Lattice, Trap>::find_groundstate(){
	for(int i = 0; i<ul_; ++i){
            s.diag(result_.begin());
            transform(output_.begin(), output_.end(), result_.begin(), output_.begin(), minus<double>());
            if (check(output_.begin(), output_.end(), result_.begin())){
                if(wi_){cerr<<it_<<" iteration(s) needed"<<endl;}
                return result_;
            } 
            else {if(wi_){++it_;}
                copy(result_.begin(), result_.end(), output_.begin() );
                s.re_init(result_, bec::Potential<T, Lattice, Trap>());
            }	
        }
    //if we arrive here, there is nothing more to do, we can only return the result.
    return result_;
}


template<class T, class MATRIX, class VECTOR, class DiagPolicy, bool Lattice, bool Trap>
 template<class In>
	bool Hamiltonian<T, MATRIX, VECTOR, DiagPolicy, Lattice, Trap>::check(In first, In last, In newvec){
		while(first != last){
			if(abs(*first)<=100.*numeric_limits<double>::epsilon()){++first, ++newvec;}
			else{return false;}
			
		}	
		//if we arrive here, everything is fine
		return true;
	}
//finished with Hamiltonian.
	
#endif //BEC_GROUNDSTATE_H
