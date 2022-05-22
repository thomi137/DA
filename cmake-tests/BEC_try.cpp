#include "../headers/BEC_Groundstate.h"
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace boost::numeric::ublas;

int main(){
	const double pi = 4*atan(1.0);
	int N = 100;
	double L = 10.;
	typedef std::vector<double> Vector;
	typedef boost::numeric::ublas::matrix<double> Matrix;
			
        Vector data(N, 1./sqrt(2.*pi) );
	double g = 0.0; 
	/*
	for(size_t i = 0; i<N; ++i){
            double xpos =L*0.5-double(i)*L/double(N);
		data[i]=1./sqrt(2*pi)*exp(-xpos*xpos/2.);
	}
	*/
	copy(data.begin(), data.end(), ostream_iterator<double>(cout, "\n") );
	cout<<"&"<<endl;
	
	Hamiltonian<Matrix, Vector, FullDiag> H(N, L, g, data, 5000, true);
	
	Vector result = H.find_groundstate();
	copy(result.begin(), result.end(), ostream_iterator<double>(cout,"\n"));
	
	return 0;
	
}
	
