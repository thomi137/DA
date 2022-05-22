#include<iostream>
#include<cmath>
#include<vector>
#include<iterator>
#include<complex>
#include<algorithm>
#include"fftw3.h"

using namespace std;

const double pi = 3.141592654;
const double twopi = 2*pi;
 

int main(){
  
  int N = 32;
  //initializing BEC ground state. 
  double theta = 0.;
  vector<complex<double> > psi(N, 1./sqrt(twopi) );
  complex<double>* it = psi.begin();
  vector<complex<double> > phi(N);
  
  copy(psi.begin(), psi.end(), ostream_iterator<complex<double> >(cout, "\n") );
  
  fftw_complex *in = new fftw_complex[psi.size()];
  fftw_complex *out = new fftw_complex[phi.size()];
  in = reinterpret_cast<fftw_complex*>(it);

  fftw_plan p;
  p=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
  
  cout<<"starting Transformation:\n";
  fftw_execute(p);
  fftw_destroy_plan(p);

  for(int i=0; i<phi.size(); ++i){
    complex<double> temp(out[i][0],out[i][1]);
    phi[i]=temp;
  }
  

  //backwards for a check:
  cout<<"Transforming back...\n";
  complex<double>* back = phi.begin();
  
  in = reinterpret_cast<fftw_complex*>(back);

  fftw_plan p1;
  p1=fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_MEASURE);

  fftw_execute(p1);
  fftw_destroy_plan(p1);

  for(int i=0; i<psi.size(); ++i){
    complex<double> temp(in[i][0]/32.,in[i][1]/32.);
    psi[i]=temp;
  }
  
  copy(psi.begin(), psi.end(), ostream_iterator<complex<double> >(cout, "\n") );

  delete[] in;
  delete[] out;
  
  return 0;
}
  
