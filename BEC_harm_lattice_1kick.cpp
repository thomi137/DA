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
#include"headers/split_step.h"
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
    
    int N=1024;
			
    const int period = 1000;   	
    const double deltat = 2.*pi/double(period);
    const double omega[4] = {1.3151861047504085, -1.1776799841887, 0.2355733213359357, 0.784513610477560};
    
    typedef complex<double> complex_type;
    typedef fftw_complex    complex_c_type;
   	
    const complex_type I(0,1);

    Vector groundstate(N, 100.);
    Cvector data(N);
    Cvector mom(N);
    Hamiltonian<double,Matrix, Vector, FullDiag, true> GP(N, L, g, groundstate, 50000, true);
    groundstate=GP.find_groundstate();
    copy(groundstate.begin(), groundstate.end(), data.begin());
    cout<<"System size"<<'\t'<<"Points"<<'\t'<<"Coupling"<<'\t'<<"deltat"<<'\t'<<"Kick strength"<<endl;
    cout<<L<<'\t'<<N<<'\t'<<g<<'\t'<<deltat<<'\t'<<K<<endl;
    
        vector<complex_type> kick(N);	
	for (size_t i = 0; i<N; ++i) {
		double xpos = L*0.5-double(i)*L/double(N);
          	kick[i]= exp(-I*K*cos(xpos));
        }
    mom=FourierTransform<fft::forward, fft::Estimate>(data);
    for(int i = 0; i<mom.size(); ++i){
	    int index = i<N/2?i+N/2:i-N/2;
	    cout<<norm(mom[index])<<endl;
    }
	    	
	transform(data.begin(), data.end(), kick.begin(), data.begin(), multiplies<complex_type>() );
	copy(data.begin(), data.end(), ostream_iterator<complex_type>(cout, "\n") );		
	
   	fftw_complex* in  = new fftw_complex[sizeof(fftw_complex)*N];
        fftw_complex* out = new fftw_complex[sizeof(fftw_complex)*N];
	
        double t=0;
        int stepcounter=0;
       /*
	for(;t<10; t+=deltat, ++stepcounter){
                
		
		data=S3(L, deltat, omega[3], g, data, in, out);
		data=S3(L, deltat, omega[2], g, data, in, out);
		data=S3(L, deltat, omega[1], g, data, in, out);
		data=S3(L, deltat, omega[0], g, data, in, out);
		data=S3(L, deltat, omega[1], g, data, in, out);
		data=S3(L, deltat, omega[2], g, data, in, out);
		data=S3(L, deltat, omega[3], g, data, in, out);
	
		double mynorm = norm_of_vec(data.begin(), data.end(), 0.);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex_type>(), sqrt(mynorm) ) );
	
               // if(!(stepcounter%period)){
                    cout<<"&"<<"\t"<<"&"<<endl;
		    bec::momentum_distribution(data.begin(), data.end(), mom.begin());
                    for(int i=0; i<data.size(); ++i)
                        cout<<norm(data[i])<<"\t"<<norm(mom[i])<<endl;
              //  }
	}*/

	return 0;
}
		
	
