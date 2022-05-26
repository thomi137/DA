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
    const double deltat = 2.*pi/double(period);
    const complex_type I(0,1);
    
    double t=0;
    int stepcounter=1;
		
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

    //finding Groundstate
    Cvector data(N);
    Cvector mom(N);
    Cvector mean(N);
      
    for(int i = 0; i<N; ++i){
        double xpos = L*0.5-double(i)*L/double(N);
        data[i]=1./sqrt(sqrt(pi) )*exp(-xpos*xpos/2.);
    }
    
    for(int i = 0; i<20000;++i){
        data=split_step::S3imag(L, deltat, g, data, bec::Potential<double, false, false>());
        mynorm = L_2_norm(data.begin(), data.end(), L);
        //cout<<mynorm<<endl;
        transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), mynorm) );
    }

    //properly normalize the whole thing, so that we may use it with the propagator.
    mynorm = L_2_norm(data.begin(), data.end(), L);
    transform(data.begin(), data.end(), data.begin(), bind2nd(divides<complex<double> >(), mynorm) );

    cout<<"System size"<<'\t'<<"Points"<<'\t'<<"Coupling"<<'\t'<<"deltat"<<'\t'<<"Kick strength"<<'\t'<<"Max time"<<'\t'<<"Period"<<endl;
    cout<<L<<'\t'<<N<<'\t'<<g<<'\t'<<deltat<<'\t'<<K<<'\t'<<maxt<<'\t'<<period<<endl<<endl;
   
       vector<complex_type> kick(N);	
	for (size_t i = 0; i<N; ++i) {
		double xpos = L*0.5-double(i)*L/double(N);
          	kick[i]= exp(-I*K*sin(xpos));
        }
	
	transform(data.begin(), data.end(), kick.begin(), data.begin(), multiplies<complex_type>() );
       
	for(;t<maxt; t+=deltat, ++stepcounter){
               
		data=split_step::S3(L, deltat, g, data, bec::Potential<double, false, false>());
		mom=FourierTransform<fft::forward, fft::Estimate>(data);
                transform(mom.begin(), mom.end(), mom.begin(), bind2nd(divides<complex<double> >(), double(period) ) );        
		transform(mean.begin(), mean.end(), mom.begin(), mean.begin(), plus<complex<double> >() );
	
                if(!(stepcounter%period) ){
                    cout<<"& "<<"period: "<<period<<endl;		    	
                    for(int i = 0; i<mean.size(); ++i){
                        int index = (i<mean.size()/2)?i+mean.size()/2:i-mean.size()/2;
                        cout<<norm(mean[index])<<endl;
                        mean[index]=0;
                    }
                }
	}
	return 0;

} 
		
	
