#include"fftw3.h"
#include<vector>
#include<iostream>
#include"headers/split_step.h"
#include"headers/bec.h"
#include<algorithm>
#include<functional>
#include<iterator>
#include<cmath>
#include<complex>
#include"headers/BEC_Groundstate.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

int main(int argc, char* argv[]){
	typedef std::vector<double> Vector;
    	typedef std::vector<complex<double> > Cvector;
    	typedef boost::numeric::ublas::matrix<double> Matrix;
    	typedef complex<double> complex_type;
    	typedef fftw_complex    complex_c_type;
    	typedef vector<complex<double> >::size_type size_t;      
    	typedef complex<double> complex_type;
    	typedef fftw_complex complex_c_type;

    	double K = 0.1;
    	double g = 0.1;
    	double maxt = 20;
    	int period = 10000;
    	const double L = 10;
    	double mynorm= 0;
    	int N=1024;
   	const complex_type I(0,1);
    
    	double t=0;
    	int stepcounter=0;
		
	//parsing command line
    	switch(argc){
		case 1:
			break;
		case 2:
			K=atof(argv[1]);
			break;
		case 3:
			K=atof(argv[1]);
			g=atof(argv[2]);
			break;
		case 4:
			K=atof(argv[1]);
			g=atof(argv[2]);
			maxt=atof(argv[3]);	
			break;
		case 5:
			K=atof(argv[1]);
			g=atof(argv[2]);
			maxt=atof(argv[3]);
			period = atoi(argv[4]);	
			break;	
			
		case 6:
			K=atof(argv[1]);
			g=atof(argv[2]);
			maxt=atof(argv[3]);
			period = atoi(argv[4]);	
			N=atoi(argv[5]);	
			break;
			
		default:
			break;
    
	}
   	double deltat = 2.*pi/double(period);
   	double alpha = double(N)/(double(N)+1.);
		
   	 //finding Groundstate
    	Cvector data(N);
    	Cvector mom(N);
    	Cvector mean(N);
      
   	 for(int i = 0; i<N; ++i){
        	double xpos = L*0.5-double(i)*L/double(N);
        	data[i]=1./sqrt(sqrt(pi) )*exp(-xpos*xpos/2.);
    	}
    
    	for(int i = 0; i<20000;++i){
        	data=split_step::S3imag(L, deltat, g, data, bec::Potential<double, true, false>());
        	mynorm = L_2_norm(data.begin(), data.end(), L);
        	//cout<<mynorm<<endl;
        	transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), mynorm) );
    	}

    	//properly normalize the whole thing, so that we may use it with the propagator.
    	mynorm = L_2_norm(data.begin(), data.end(), L);
    	transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), mynorm) );

    	cout<<"System size"<<'\t'<<"Points"<<'\t'<<"Coupling"<<'\t'<<"deltat"<<'\t'<<"Kick strength"<<'\t'<<"Max time"<<'\t'<<"Period"<<endl;
    	cout<<L<<'\t'<<N<<'\t'<<g<<'\t'<<deltat<<'\t'<<K<<'\t'<<maxt<<'\t'<<period<<endl<<endl;
   
  	//Good. We're at the ground state of the one-dimensional condensate. We can now proceed with the kick...
    	
	
       	vector<complex_type> kick(N);	
	for (size_t i = 0; i<N; ++i) {
		double xpos = L*0.5-double(i)*L/double(N);
          	kick[i]= exp(-I*K*sin(xpos));
        }
					
    	//So much for the kick. Now, if the kick is strong enough, we expect it to yield a rapid phase decoherence.
    	//However, when integrating the GP equation over large times we necessarily catch round-offs, hence some noise.
    	//It is therefore much more practical to look only at some of the timesteps. So we first catch a reasonable mean and 
    	//proceed by taking the exponentially weighted moving average, printing it out after some times.
    	//Details (and lembas...) can be found at http://lorien.ncl.ac.uk/ming/filter/filewma.htm.
	
	while( stepcounter++ < (period - 1)){
		data=split_step::S3(L, deltat, g, data, bec::Potential<double, true, false>());
		mom=FourierTransform<fft::forward, fft::Estimate>(data);
		transform(mom.begin(), mom.end(), mom.begin(), bind2nd(divides<complex<double> >(), double(period)) );
		transform(mean.begin(), mean.end(), mom.begin(), mean.begin(), plus<complex<double> >() );
	}
    
	//We have now calculated sort of a moving average over *period* distributions. We will now proceed weighting it with
	//exponential weights, the most recent distribution being given the biggest weight. Note that this *period*
	//distributions mark one unit of time. This of course also means that the first few distributions gain equal weigth, but 
	//this should be ok, since here the most physics should take place...
	
	++stepcounter;
	for(;t<maxt; t+=deltat){
               	
		data=split_step::S3(L, deltat, g, data, bec::Potential<double, true, false>());
		mom=FourierTransform<fft::forward, fft::Estimate>(data);
		
		if(stepcounter==(period*50+1)){
			cerr<<"Kick at time: "<<t<<endl;
			transform(data.begin(), data.end(), kick.begin(), data.begin(), multiplies<complex_type>() );
		}
		
		//we take the weights according to X(n)=alpha*X(n-1)+(1-alpha)*x and take this to be 
		//the filtered value of the distribution at time t.
		
		transform(mean.begin(), mean.end(), mean.begin(), bind2nd(multiplies<complex<double> >(), alpha) );
                transform(mom.begin(), mom.end(), mom.begin(), bind2nd(multiplies<complex<double> >(), (1.-alpha) ) );       
		transform(mean.begin(), mean.end(), mom.begin(), mean.begin(), plus<complex<double> >() );
	
		//and sometimes, we even print the whole stuff out...
                if(!(stepcounter%period) ){
                    cout<<"& "<<"Step: "<<t<<endl;		    	
                    for(int i = 0; i<mean.size(); ++i){
                        int index = (i<mean.size()/2)?i+mean.size()/2:i-mean.size()/2;
                        cout<<norm(mean[index])<<endl;
                    }
                }
		++stepcounter;
	}
	return 0;

} 
		
	
