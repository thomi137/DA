#include"headers/fft.h"
#include<complex>
#include<vector>
#include<algorithm>
#include<iostream>
#include<iterator>
#include<cmath>

using namespace std;

int main(){
	double pi = 4*atan(1.0);
	int N = 500;
	typedef complex<double> cdouble;
	vector<cdouble> x(100);
	vector<cdouble> p(100);

	for(int i = 0; i < N; ++i){
		double xpos = 5. - i*10./double(N);
		x[i]=exp(-xpos*xpos);
}

p = FourierTransform<fft::forward, fft::Estimate>(x);

	for(int i = 0; i < N; ++i){
		double xpos = 5.-i*10./double(N);    
		int index = (i<N/2.)?i+N/2:i-N/2;
		cout<<xpos<<'\t'<<abs(x[i])<<'\t'<<sqrt(2*pi)*abs(p[index])<<endl;
}


return 0;
}
