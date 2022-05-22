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

inline vector<complex<double> > S3 (const double deltat, const double omega, const double g, vector<complex<double>
		> data, const vector<complex<double> >& potential, fftw_complex* in, fftw_complex* out)
		{
			complex<double>I(0,1);
			unsigned int N = data.size();
			vector<complex<double> > kvec(N);
			vector<complex<double> > kfac(N);
			
			for (unsigned int i = 0; i<N; ++i) {
				kvec[i]= (double(i)<double(N)/2.) ? double(i): -double(N-i);
          			kfac[i]= exp(-1.*deltat*omega/2.*kvec[i]*kvec[i]/2.);
			}
						
			in=reinterpret_cast<fftw_complex*>(&data.front());
                	fftw_plan fourier1 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(fourier1);
			fftw_destroy_plan(fourier1);
                	for(unsigned int i = 0; i<N; ++i){
	                	double real_part = out[i][0]; double imag_part = out[i][1];
				complex<double> fac = kfac[i];
                        	out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
				out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
			}
                	fftw_plan backfourier1 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	        	fftw_execute(backfourier1);
			fftw_destroy_plan(backfourier1);	
			
			for(unsigned int i = 0; i<N; ++i) 
                		  data[i]*=exp(-1*g*deltat*norm(data[i])); 
			
			transform(data.begin(), data.end(), potential.begin(), data.begin(), multiplies<complex<double> >() );
				
			in=reinterpret_cast<fftw_complex*>(&data.front());
			fftw_plan fourier3 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(fourier3);
			fftw_destroy_plan(fourier3);
                		for(unsigned int i = 0; i<N; ++i){
	                	double real_part = out[i][0]; double imag_part = out[i][1];
				complex<double> fac = kfac[i];
                        	out[i][0] = (real(fac)*real_part-imag(fac)*imag_part)/double(N);
				out[i][1] = (real(fac)*imag_part+imag(fac)*real_part)/double(N);
			}
                	fftw_plan backfourier3 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	        	fftw_execute(backfourier3);
			fftw_destroy_plan(backfourier3);
			return data;				
		}

int main(){
	const double pi = 4.*atan(1.);
        const double g = 0.5;
	const double L = 10.;
	const double q = 40.*pi/L;
	double newnorm = 0.;
	double d = 0.1;
	int N=1024;
			
	const int period = 10000;	
	const double deltat = 2.*pi/double(period);
	const double omega[4] = {1.3151861047504085, -1.1776799841887, 0.2355733213359357, 0.784513610477560};
	
	typedef complex<double> complex_type;
	typedef fftw_complex    complex_c_type;
	typedef vector<complex<double> >::size_type size_t;
	
	const complex_type I(0,1);

        vector<complex_type> data(N, 100.);
	/*
	for(size_t i = 0; i<N; ++i){
            double xpos =20.*0.5-double(i)*20./double(N);
		data[i]=1./sqrt(2*pi)*exp(-xpos*xpos/2.);
	}
		*/
	vector<complex<double> > potential(N);
		for(unsigned int i = 0; i<N; ++i){
			double xpos =L*0.5-double(i)*L/double(N);
			double sinqx = sin(q*xpos);
			potential[i]=exp(-1.*deltat*((xpos*xpos/2.)+q*q*0.5*sinqx*sinqx));
		}

   	fftw_complex* in  = new fftw_complex[sizeof(fftw_complex)*N];
	fftw_complex* out = new fftw_complex[sizeof(fftw_complex)*N];
	
        double t=0;
        int stepcounter=0;
	for(;t<50; t+=deltat, ++stepcounter){
	
		
		data=S3(deltat, 1., g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		
		/*
		data=S3(deltat, omega[3], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[2], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[1], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[0], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[1], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[2], g, data, potential, in, out);
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		data=S3(deltat, omega[3], g, data, potential, in, out);		
		newnorm = L_2_norm(data.begin(), data.end(), L);
		transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), newnorm)	);
		*/
		
		
	
	}
	
	
	for(int i = 0; i<N;++i)
			cout<<norm(data[i])<<endl;
		cout<<"&"<<endl;	
	
	return 0;
} 
		
	
