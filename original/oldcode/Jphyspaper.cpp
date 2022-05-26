#include"bec.h"
#include<vector>
#include<iterator>
#include<algorithm>
#include<iostream>
#include<cmath>

using namespace std;

int main(){
	
	const complex<double> I(0,1);
  	double g = 1.;
  	const double deltat = 0.0025;
  	int N=1024;
  
 	//k (momentum space)
  	vector<complex<double> > k(N);
  	for(unsigned int i = 0; i < k.size(); ++i)
      		k[i] = (i<N/2.) ? 2*M_PI/260.*double(i):-2*M_PI/260.*double(N-i);
  
  	vector<complex<double> > data(N);
  	vector<complex<double> > data_copy(data.size());
  	vector<complex<double> >::pointer data_pointer = &data.front();
  	vector<complex<double> >::pointer data_copy_pointer = &data_copy.front();
  
	//Soliton solution.
	for(int i = 0; i < N; ++i){
      		double xposition =  -20.+260./double(N)*double(i);
      		data[i]=exp(-I*(2.*xposition+M_PI/2.))*1./cosh(xposition);
  	}

	//prints the soliton solution squared.
	for(int i = 0; i<N; ++i)
    	cout<<norm(data[i])<<endl;
  
	//fftw related things	.	
	//some arrays to store data in.
	fftw_complex* in  = new fftw_complex[sizeof(fftw_complex)*N];
	fftw_complex* out = new fftw_complex[sizeof(fftw_complex)*N];
	
	//propagation
	for(double t = 0.;t<=50.;t+=deltat){
    
   	 	//get a copy of data.
    		copy(data.begin(), data.end(), data_copy.begin());

    		//point the input array to data_copy.
      		in  = reinterpret_cast<fftw_complex*>(data_copy_pointer);
    		fftw_plan nonlinear_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    		//transform the amplitude to k-space
    		fftw_execute(nonlinear_plan);
    		fftw_destroy_plan(nonlinear_plan);

    		//getting exp(dx^2...)*psi
    		for(int i = 0; i < N; ++i){
    			double real_part = out[i][0];
			double imag_part = out[i][1];
        		complex<double> exp_fac = exp(-I*deltat*(cos(k[i])-1));
       			out[i][0]=(real(exp_fac)*real_part-imag(exp_fac)*imag_part)/double(N);
      			out[i][1]=(real(exp_fac)*imag_part+imag(exp_fac)*real_part)/double(N);
		}

   		//transforming back.
    		fftw_plan nonlinear_backtransform  = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);	
    		fftw_execute(nonlinear_backtransform);
    		fftw_destroy_plan(nonlinear_backtransform);
 
		//nonlinear local interaction (in real space...)
    		for(unsigned int i = 0; i<data.size(); ++i){
  			complex<double> temp = conj(data_copy[i])*data_copy[i]; 
      			data[i]*=exp(I*g*deltat*temp);
    		}

    
   	 	//point the second input array to data.
    		in = reinterpret_cast<fftw_complex*>(data_pointer);
    		fftw_plan kinetic_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		
    		fftw_execute(kinetic_plan);
    		fftw_destroy_plan(kinetic_plan);

		//kinetic term
    		for(int i=0; i<N;++i){
     			double real_part = out[i][0];
      			double imag_part = out[i][1];
      			complex<double> fac = exp(-I*deltat*2*(cos(k[i])-1));
      			out[i][0]=(real_part*real(fac)-imag_part*imag(fac))/double(N);
      			out[i][1]=(real_part*imag(fac)+imag_part*real(fac))/double(N);
    		}

   		fftw_plan kinetic_backtransform = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    		fftw_execute(kinetic_backtransform);
    		fftw_destroy_plan(kinetic_backtransform);
 	}
	


  	cout<<"&"<<endl;
  	for(int i=0;i<N;++i)
      		cout<<norm(data[i])<<endl;

  //cleaning up...
  //fftw_destroy_plan(new_plan);
  //fftw_destroy_plan(backplan);
  fftw_cleanup();
  //*/ 
  return 0;
  /*
    for(int i = 0; i<N; ++i){
        cout<<in[i][0]<<","<<in[i][1]<<endl;
    }
*/
		
		
/*  cout<<endl;

    for(int i = 0; i<N; ++i){
        cout<<data[i]<<endl;
    }

    cout<<endl;
 */   
    
}
