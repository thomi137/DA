#include"fftw3.h"
#include<vector>
#include<iostream>
#include"headers/bec.h"
#include<algorithm>
#include<functional>
#include<iterator>
#include<cmath>
#include<complex>

using namespace std;

int main(){
	const double pi = 4.*atan(1.);
	const double K = 0.8;
        const double g = 2.1;
	int N=1024;
			
	const int period = 1000;	
	const double deltat = 2.*pi/double(period);
	
	typedef complex<double> complex_type;
	typedef fftw_complex    complex_c_type;
	typedef vector<complex<double> >::size_type size_t;
	
	const complex_type I(0,1);

        vector<complex_type> data(N, 1./sqrt(2.*pi));	
	vector<complex_type> kfac(N);
	vector<complex_type> kick(N);
		
	for (size_t i = 0; i<N; ++i) {
          kfac[i]= (double(i)<double(N)/2.)? double(i): -1.*double(N-i);
          kfac[i]= exp(-I*deltat/2.*kfac[i]*kfac[i]/2.);
          kick[i]= exp(-I*K*cos(2.*pi/double(N)*double(i)));
        }
	
	complex_c_type* in  = new complex_c_type[sizeof(complex_c_type)*N];
	complex_c_type* out = new complex_c_type[sizeof(complex_c_type)*N];
        in = reinterpret_cast<complex_c_type*>(&data.front());
	
        double t=0;
        int stepcounter=0;
	for(;t<1000; t+=deltat, ++stepcounter){

		if(!(stepcounter%period))
                 transform(data.begin(), data.end(), kick.begin(), data.begin(), multiplies<complex_type>() );
		
                fftw_plan fourier1 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(fourier1);
                for(int i = 0; i<N; ++i){
	                double real_part = out[i][0]; double imag_part = out[i][1];
			complex_type fac = kfac[i];
                        out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
			out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
		}
                fftw_plan backfourier1 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	        fftw_execute(backfourier1);	

                for(size_t i = 0; i<data.size(); ++i) 
                  data[i]*=exp(-I*g*deltat*norm(data[i])); //something wrong here?
                
                fftw_plan fourier2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(fourier2);
                for(int i = 0; i<N; ++i){
			double real_part = out[i][0]; double imag_part = out[i][1];
			complex_type fac = kfac[i];
                        out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
			out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
		}
                fftw_plan backfourier2 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(backfourier2);

		double mynorm = norm_of_vec(data.begin(), data.end(), 0.);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex_type>(), sqrt(mynorm) ) );
		
		if(!(stepcounter%period))
                  cout<<" " << real(bec::energy(data, g))<<endl;                 
	}
	return 0;
} 
		
	
