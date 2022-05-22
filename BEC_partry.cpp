#include<iostream>
#include<cmath>
#include<vector>
#include<iterator>
#include<algorithm>
#include<limits>
#include<functional>
#include<boost/random/lagged_fibonacci.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;

int main(){
	typedef std::vector<double> Vector;
	typedef boost::numeric::ublas::matrix<double> Matrix;
	
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


	
	return 0;
	
}
	
