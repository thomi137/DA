#include"bec.h"
#include<vector>
#include<iterator>
#include<algorithm>
#include<functional>
#include<iostream>
#include<cmath>


using namespace std;



int main(){

  const double pi = M_PI;
  const complex<double> I(0,1);
  double t = 0.;
  double g = 0.0;
  double timesteps = 100.;
  int stepcounter = 0;
  const double deltat = 2*pi/timesteps;
  int N=128;
  double K = 0.8;
  bool output=false;

  //kick vector
  vector<complex<double> > kick(N);
  for(unsigned int i = 0; i < kick.size(); ++i)
    kick[i]=exp(-I*deltat*K*cos(2*pi/N*i) );
  
  //k vector
  vector<complex<double> > k(N);
  for(unsigned int i = 0; i < k.size(); ++i)
    k[i] = (i<N/2.)? double(N)*double(i) : -1.*double(N)*double(N-i);

  //Initialize BEC Groundstate.
  vector<complex<double> > data(N, 1./sqrt(2.*pi));
  vector<complex<double> > data_copy(data);

  //the other stuff...
  //exponential factor in fourier space:
  vector<complex<double> >exp_k_four_space(N);
  for(unsigned int i = 0; i < exp_k_four_space.size(); ++i)
    exp_k_four_space[i] = exp(I*deltat*(cos(k[i])-1));

  //fftw related things.
  fftw_complex *in  = new fftw_complex[N];
  fftw_complex *out = new fftw_complex[N];
  
  while(t<=10000.){
  
    //get a copy   
    //it = data_copy.begin();
   // copy(data.begin(), data.end(), data_copy.begin());

    //point the input array to data_copy.
    //in = reinterpret_cast<fftw_complex*>(point);

    //Split Operator technique: multiply by the splitted operator:
    
    if (!(stepcounter%100)){
      	for(unsigned int i = 0; i < data.size(); ++i)
		data[i]*=kick[i];
      	output=true;
    }
  /* 
    fftw_plan p;
    p=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

     //getting the exponent.
     for(int i = 0; i < N; ++i){
       complex<double> exp_fac = exp_k_four_space[i];
       out[i][0]=(real(exp_fac)*out[i][0]-imag(exp_fac)*out[i][1])/double(N);
       out[i][1]=(real(exp_fac)*out[i][1]+imag(exp_fac)*out[i][0])/double(N);
     }

     //the above is now written to data_copy
     fftw_plan bp;
     bp=fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
     fftw_execute(bp);
   
     //nonlinear interaction:
     for(int i = 0; i<data.size(); ++i){
       //       cout<<data[i]<<endl;
       data[i]*=exp(-I*g*deltat*temp);
       //       cout<<data[i]<<endl;
     }
*/
     //preparing for derivative...
     in=reinterpret_cast<fftw_complex*>(&data.front());
     fftw_plan kinetic_transform=fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
     fftw_execute(kinetic_transform);
     fftw_destroy_plan(kinetic_transform);
  
    //kinetic term
    for(int i=0; i<N;++i){
      double real_part = out[i][0];
      double imag_part = out[i][1];
      complex<double> fac = exp_k_four_space[i];
      out[i][0]=(real_part*real(fac)-imag_part*imag(fac))/double(N);
      out[i][1]=(real_part*imag(fac)+imag_part*real(fac))/double(N);
    }
   
    fftw_plan finish = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(finish);
    fftw_destroy_plan(finish);

    if(output){
   	for(unsigned int i = 0; i<data.size();++i)
		cout<<norm(data[i])<<endl;
	cout<<"&"<<endl;
    }
   
    t+=deltat;
    stepcounter++;
    output=false;

  }
   
  //cleaning up...
  //fftw_destroy_plan(bp);
  fftw_cleanup();
  //*/ 
  return 0;
}
