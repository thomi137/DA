#include"fftw3.h"
#include<vector>
#include<iostream>
#include"headers/bec.h"
#include<algorithm>
#include<functional>
#include<iterator>
#include<cmath>
#include<complex>
#include"headers/BEC_Groundstate.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

int main(int argc, char* argv[]){
    typedef std::vector<double> Vector;
    typedef std::vector<complex<double> > Cvector;
    typedef boost::numeric::ublas::matrix<double> Matrix;
    typedef complex<double> complex_type;
    typedef fftw_complex    complex_c_type;
    typedef vector<complex<double> >::size_type size_t;
    
    const double pi = 4.*atan(1.);
     double K = 0.1;
     double g = 0.1;
    const double L = 10;
    		
    int N=1024;
			
    switch(argc){
		case 1:
			break;
		case 2:
			K=atof(argv[1]);
			break;
		case 3:
			K=atof(argv[1]);
			g=atof(argv[2]);
			break;
		default:
			break;
    
		}
    
	vector<complex_type> kick(N);	
	for (size_t i = 0; i<N; ++i) {
		double xpos = L*0.5-double(i)*L/double(N);
          	kick[i]= exp(-I*K*sin(2*pi/L*xpos));
        }	

    const int period = 100;   	
    const double deltat = 2.*pi/double(period);
  
    typedef complex<double> complex_type;
    typedef fftw_complex    complex_c_type;
   	

    Vector groundstate(N, 100.);
    Hamiltonian<double, Matrix, Vector, FullDiag, true> GP(N, L, g, groundstate, 50000, true);
    groundstate=GP.find_groundstate();
    vector<complex_type> mom(N);
  vector<complex<double> > data(groundstate.size() );
  copy(groundstate.begin(), groundstate.end(), data.begin() );
    copy(data.begin(), data.end(), ostream_iterator<complex<double> >(cout, "\n") );
        
  
    transform(data.begin(), data.end(), kick.begin(), data.begin(), multiplies<complex<double> >() );

    cout<<endl;
    
    bec::momentum_distribution(data.begin(), data.end(), mom.begin());
    for (int i=0; i<mom.size();++i)
	cout<<norm(mom[i])<<endl;

	return 0;
} 
		
	
