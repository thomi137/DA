#include<complex>
#include"fftw3.h"
#include<vector>
#include<iterator>
#include<algorithm>
#include<functional>
#include<iostream>
#include<cmath>

using namespace std;


int main(){

  const double pi = M_PI;
  const double twopi = 2.*M_PI;
  const complex<double> I(0,1);
  double t = 0;
  double g = 1.5;
  double timesteps = 10.;
  bool output = false;
  int stepcounter = 0;
  const double deltat = twopi/timesteps;
  int N=1024;
  double K = 0.8;
  //kick vector
  vector<complex<double> > kick(N);
  for(int i = 0; i < kick.size(); ++i)
    kick[i]=exp(-I*deltat*K*cos(twopi/N*i) );
  //Initialize BEC Groundstate.
  vector<complex<double> > data(N, 1./sqrt(twopi));
 // copy(data.begin(), data.end(), ostream_iterator<complex<double> >(cout,"\n") );
 // cout<<endl;


  //the other stuff...
  vector<complex<double> >::pointer point = &data.front();
  vector<complex<double> >result(N);
  vector<complex<double> >::pointer point1 = &result.front();
  fftw_complex *in  = new fftw_complex[sizeof(fftw_complex)*N];
  fftw_complex *out = new fftw_complex[sizeof(fftw_complex)*N];
  complex<double>* helper1 = new complex<double>[sizeof(complex<double>) * N];
  
 
  while(t<=1000){
    //Split Operator technique: multiply by the splitted operator:
    if (!(stepcounter%10)){//if so, then kick!
     output=true;
      for(int i = 0; i < data.size(); ++i)
	data[i]*=kick[i];
    }
    //then do the rest...
    
    //Debugging Code
    // copy(data.begin(), data.end(), ostream_iterator<complex<double> >(cout,"\n") );
    
    in = reinterpret_cast<fftw_complex* >(point);
   
    fftw_plan p;
    p=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    for(int i = 0; i<N;++i){
      double real_part = out[i][0];
      double imag_part = out[i][1];
      complex<double> fac = exp(-1*I*deltat*(cos((i<N/2.)?double(N)*i:-1*double(N)*double(N-i))-1.));
      out[i][0]=real_part*real(fac)+imag_part*imag(fac);
      out[i][1]=real_part*imag(fac)+imag_part*real(fac);
    }
    
    fftw_plan bp;
    bp=fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(bp);
    fftw_destroy_plan(bp);
    helper1=reinterpret_cast<complex<double>* >(in);
    copy(helper1, helper1+N-1, data.begin() );
    for(int i = 0; i<data.size(); ++i)
      data[i]/=N;

    if(output){
for (int i = 0; i<N;++i)
		cout<<norm(data[i])<<endl;
cout<<"&"<<endl;
}
    cout<<endl;
    t+=deltat;
    stepcounter++;
    
  }
//*/
  return 0;
}
