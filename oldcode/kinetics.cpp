#include<complex>
#include"fftw3.h"
#include<vector>
#include<iostream>

using namespace std;

int main(){
	int N = 4;
	double deltat = 0.0125;
	const complex<double> I(0,1);
	vector<complex<double> >vec(N,1./sqrt(2.*M_PI) );		
	fftw_complex* in = new fftw_complex[N];
	in = reinterpret_cast<fftw_complex*>(&vec.front() );
	fftw_plan p = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	
	copy(vec.begin(), vec.end(), ostream_iterator<complex<double> >(cout, "\n") );
	cout<<endl;

	for(int i = 0; i<N; ++i){
		double re=in[i][0];
		double im=in[i][1];
		complex<double> kvec = (i<N/2.) ? I*double(i): -I*double(N-i);
		complex<double> fac = exp(-I*deltat*kvec*kvec);
		in[i][0] = (real(fac)*re-imag(fac)*im)/double(N);
		in[i][1] = (real(fac)*im+imag(fac)*re)/double(N);
		}
	
	fftw_execute(p);
	fftw_destroy_plan(p);
		
	copy(vec.begin(), vec.end(), ostream_iterator<complex<double> >(cout, "\n") );
	
	return 0;
	}
