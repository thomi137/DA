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
  int stepcounter = 0;
  const double deltat = twopi/timesteps;
  int N=1024;
  double K = 0.8;
  bool output;

  //kick vector
  vector<complex<double> > kick(N);
  vector<complex<double> > k(N);
  for(int i = 0; i < kick.size(); ++i)
    kick[i]=exp(-I*deltat*K*cos(twopi/N*i) );
  for(int i = 0; i < k.size(); ++i)
    k[i] = (i< N/2.)? I*double(i) : -I*double(N-i);

  //Initialize BEC Groundstate.
  vector<complex<double> > data(N, 1./sqrt(twopi));
  for (int i = 0; i < data.size(); ++i)
    cout<<2*pi/N*i<<'\t'<<abs(data[i])<<endl;
  cout<<"&"<<endl;
  


  //the other stuff...
  vector<complex<double> >::pointer kpoint = &k.front();
  vector<complex<double> >::pointer point  = &data.front();
  vector<complex<double> >result(N);
  vector<complex<double> >exp_k_real_space(N);
  fftw_complex *in  = new fftw_complex[sizeof(fftw_complex)*N];
  fftw_complex *out = new fftw_complex[sizeof(fftw_complex)*N];
  complex<double>* helper1 = new complex<double>[sizeof(complex<double>) * N];
  
  in = reinterpret_cast<fftw_complex* >(kpoint);

  fftw_plan bp;
  bp=fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(bp);
 
  helper1=reinterpret_cast<complex<double>* >(out);
  copy(helper1, helper1+N-1, exp_k_real_space.begin());
  copy(exp_k_real_space.begin(), exp_k_real_space.end(), ostream_iterator<complex<double> > (cout, "\n"));
    
/*
   while(t<=1000){
     //Split Operator technique: multiply by the splitted operator:
     if (!(stepcounter%10)){//if so, then kick!
     // cout<<"Now Kicking..."<<endl;
      for(int i = 0; i < data.size(); ++i)
	data[i]*=kick[i];
      output=true;
     }
    //then do the rest...
      


    fftw_plan p;
    p=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    
    for(int i=0; i<N;++i){
      double real_part = out[i][0];
      double imag_part = out[i][1];
      complex<double> fac = exp(I*deltat*k[i]*k[i]/2.);
      out[i][0]=real_part*real(fac)-imag_part*imag(fac);
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
    for (int i = 0; i < data.size(); ++i)
   cout<<2*pi/N*i<<'\t'<<abs(data[i])<<endl; 

    cout<<"&"<<endl;
    }
    t+=deltat;
    stepcounter++;
    output=false;
    }
//*/
  return 0;
}
