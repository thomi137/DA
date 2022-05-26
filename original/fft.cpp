#include<complex>
#include"headers/fft.h"
#include<iostream>
#include<vector>
#include<functional>

using namespace std;

int main(){
	vector<complex<double> > vec(4);
    vec[1]=0.5;
    vec[2]=0.5;
	vector<complex<double> >ft = FourierTransform<fft::forward,
			fft::Estimate>(vec);

	copy(ft.begin(), ft.end(), ostream_iterator<complex<double> >(cout, "\n") );
        
        vector<complex<double> >backft = FourierTransform<fft::backward, fft::Estimate>(ft);
        copy(backft.begin(), backft.end(), ostream_iterator<complex<double> >(cout, "\n") );
        ft = FourierTransform<fft::forward,
            fft::Estimate>(backft);

        
  copy(ft.begin(), ft.end(), ostream_iterator<complex<double> >(cout, "\n") );
    }
