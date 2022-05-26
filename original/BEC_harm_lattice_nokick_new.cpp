#include<vector>
#include<iostream>
//#include"bec.h"
#include"headers/split_step.h"
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

int main(){
    typedef std::vector<double> Vector;
    typedef std::vector<complex<double> > Cvector;
    typedef boost::numeric::ublas::matrix<double> Matrix;
    typedef complex<double> complex_type;
    typedef fftw_complex    complex_c_type;
    typedef vector<complex<double> >::size_type size_t;
    
    const double pi = 4.*atan(1.);
    const double K = 0.8;
    const double g = 0.0;
    const double L = 15;
    
    int N=1024;
			
    const int period = 10000;   	
    const double deltat = 2.*pi/double(period);
     
    typedef complex<double> complex_type;
    typedef fftw_complex    complex_c_type;
   	
    const complex_type I(0,1);

    Vector groundstate(N, 100.);
    Cvector data(N);
    Hamiltonian<Matrix, Vector, FullDiag> GP(N, L, g, groundstate, 50000, true);
    groundstate=GP.find_groundstate();
    copy(groundstate.begin(), groundstate.end(), data.begin());
    cout<<"System size"<<'\t'<<"Points"<<'\t'<<"Coupling"<<'\t'<<"deltat"<<endl;
    cout<<L<<'\t'<<N<<'\t'<<g<<'\t'<<deltat<<endl;
    
    Cvector mom(N);
        mom= FourierTransform<fft::forward,fft::Estimate>(data);
    for(int i = 0 ; i < mom.size(); ++i)
        cout<<norm(mom[i])<<endl;
   
    /*
    double t=0;
        int stepcounter=0;
	for(;t<100; t+=deltat, ++stepcounter){
		

            data=split_step::S3(L, deltat,g, data);

	
                if(!(stepcounter%period)){
                    for(int i=0; i<data.size(); ++i)
                        cout<<t<<'\t'<<abs(data[i])<<endl;
                    cout<<"&"<<endl;
                }
	}
*/
	return 0;
} 
		
	
