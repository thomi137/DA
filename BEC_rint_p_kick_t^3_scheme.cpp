#include<iostream>
#include"headers/bec.h"
#include<algorithm>
#include<functional>
#include<iterator>
#include<cmath>

using namespace std;

int main(){
	const double pi = 4.*atan(1.);
	const double K = 0.8;
        const double g = 0.;
	int N=1024;
			
	const int period = 100;	
	int stepcounter = 0;	
	const double deltat = 2.*pi/double(period);
	double t = 0.;
	
	typedef complex<double> complex_type;
	typedef fftw_complex complex_c_type;
	typedef vector<complex<double> >::size_type size_t;
	
	const complex_type I(0,1);
	bool output = false;
	
	//BEC groundstate wavefunction.
	vector<complex_type> data(N, 1./sqrt(2.*pi));
      		
	complex_c_type* in = new complex_c_type[sizeof(complex_c_type)*N];
	complex_c_type* out = new complex_c_type[sizeof(complex_c_type)*N];
	
	vector<complex_type> k(N);
	vector<complex_type> kick(N);
	vector<complex_type> data_copy(N);
		
	for (size_t i = 0; i<k.size(); ++i)
		k[i]= (double(i)<double(N)/2.)? double(i): -1.*double(N-i);
        for (size_t i = 0; i<kick.size(); ++i)
		kick[i] = -I*K*cos(2.*pi/double(N)*double(i) );
	
	
	for(;t<=1000; t+=deltat, ++stepcounter){

		copy(data.begin(), data.end(), data_copy.begin());
		
		//kick
		if(!(stepcounter%period)){
			output = true;
			transform(data.begin(), data.end(), temp.begin(), data.begin(), plus<complex_type>() );
		}
   
		in = reinterpret_cast<complex_c_type*>(&data.front());
		fftw_plan kin = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(kin);
                
		for(int i = 0; i<N; ++i){
			double real_part = out[i][0]; double imag_part = out[i][1];
			complex_type fac = exp(-I*deltat*k[i]*k[i]/4.);
                        out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
			out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
		}
		
		fftw_plan finish1 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(finish1);
		
		fftw_destroy_plan(kin); fftw_destroy_plan(finish1);	
				
		in = reinterpret_cast<complex_c_type*>(&data_copy.front());
		fftw_plan nolin = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(nolin);
                
		for(int i = 0; i<N; ++i){
			double real_part = out[i][0]; double imag_part = out[i][1];
			complex_type fac = exp(-I*deltat*k[i]*k[i]/2.);
                        out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
			out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
		}
		
		fftw_plan finish2 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(finish2);
		
		fftw_destroy_plan(nolin); fftw_destroy_plan(finish2);	
				
		for(size_t i = 0; i<data.size(); ++i)
			data[i]*=exp(-I*g*deltat*norm(data_copy[i]) );
				
		//for the last time...
		copy(data.begin(), data.end(), data_copy.begin());
		
		in = reinterpret_cast<complex_c_type*>(&data.front());
		fftw_plan kin2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(kin2);
                
		for(int i = 0; i<N; ++i){
			double real_part = out[i][0]; double imag_part = out[i][1];
			complex_type fac = exp(-I*deltat*k[i]*k[i]/4.);
                        out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
			out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
		}
		
		fftw_plan finish3 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(finish3);
		
		fftw_destroy_plan(kin2); fftw_destroy_plan(finish3);	
				
		//double mynorm = norm_of_vec(data.begin(), data.end(), 0.);
                //transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex_type>(), sqrt(2.*mynorm) ) );
					   
                if(output){
			
                	cout<<real(bec::energy(data, g))/(2.*M_PI)<<endl;  

                    //cout<</*real(bec::energy(data, g))/(2.*M_PI)<<'\t'<<*/<<endl;		
		output=false;
		}
               
	}
	//*/
	return 0;
} 
		
	
