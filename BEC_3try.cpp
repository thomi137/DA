#include<iostream>
#include<cmath>
#include<vector>
#include<iterator>
#include<algorithm>
#include<limits>
#include<functional>
#include<ietl/lanczos.h>
#include<ietl/interface/ublas.h>
#include<ietl/vectorspace.h>
#include<boost/random/lagged_fibonacci.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace boost::numeric::ublas;

int main(){
	typedef boost::lagged_fibonacci607 Gen;
	typedef std::vector<double> Vector;
	typedef boost::numeric::ublas::matrix<double> Matrix;
	typedef ietl::vectorspace<boost::numeric::ublas::vector<double> > Vecspace;
	
	const double pi = 4*atan(1.0);
	int N = 1024;
	double L = 30.;
	double g = 0.0; 
	double q_=L/5.;
	double deltax = L/double(N);
		
	Matrix Rmat_(N,N);
	for(int i = 0; i<Rmat_.size1(); ++i){
		double xpos =L*0.5-double(i)*L/double(N);
        	double sinx =sin(q_*xpos);
		for(int j = 0; j<Rmat_.size2(); ++j){ 
			//if you eliminate the sine terms, then this is just a normal harmonic oscillator.
			Rmat_(i,j)=((i==j)?1./(deltax*deltax)+xpos*xpos*0.5+q_*q_*0.5*sinx*sinx
					: ((i==j-1)||(i==j+1))? -0.5/(deltax*deltax):0);
		}
	}
	Gen mygen;
	Vecspace vec(N);
	Vector eigen;
	std::vector<boost::numeric::ublas::vector<double> > eigenv;
	Vector err;
	std::vector<int> mult;
	
	ietl::lanczos<Matrix, Vecspace> lanczos(Rmat_, vec);
	//                                           max_iter, n_lowest, abs_tol, rel_tol (or vice versa...)
	ietl::lanczos_iteration_nlowest<double> iter(50000, 5, 500.*(numeric_limits<double>::epsilon(),
				std::pow(numeric_limits<double>::epsilon(), 2./3.) ) );
	
	try{
		lanczos.calculate_eigenvalues(iter, mygen);
	   	eigen=lanczos.eigenvalues();
	  	err=lanczos.errors();
	  	mult = lanczos.multiplicities();
  	 	}
	 	catch(std::runtime_error& e){
			 cout<<e.what()<<endl;

	 	}/*
		copy(eigen.begin(), eigen.end(), ostream_iterator<double>(cout,"\n") );
		cout<<endl;
		copy(mult.begin(), mult.end(), ostream_iterator<int>(cout,"\n") );
		cout<<endl;
		copy(err.begin(), err.end(), ostream_iterator<double>(cout,"\n") );
		cout<<endl;
*/
		Vector::iterator start = eigen.begin();
	 	Vector::iterator end = eigen.begin()+5;
		ietl::Info<double> info;

		try{
		    	lanczos.eigenvectors(start, end, back_inserter(eigenv), info, mygen);
    		}
    		catch(std::runtime_error& e){
			cout<<e.what()<<endl;
		}
		copy(eigenv.begin()->begin(), eigenv.begin()->end(), ostream_iterator<double>(cout,"\n")
				);
	
	return 0;
	
}
	
