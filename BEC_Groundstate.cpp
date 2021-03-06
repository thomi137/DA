#include <iostream>
#include <vector> 
#include <algorithm>
#include <functional>
#include <cmath>
#include <climits>
#include <iterator>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;


extern "C" void dsyev_(char* jobz, char* uplo,int* n,double* a,int* lda,double* w ,double*
		work, int* lwork, int* info);

typedef complex<double> cmpl;
typedef boost::numeric::ublas::vector<double> Vector;
typedef std::vector<double> StlVec;

  	    
int main (){
	const double pi = 4*atan(1.0);
	double L = 10.;
	int N = 1024;
	const double q = pi*40./L;
	double deltax = L/double(N);
	double g=0;

	std::vector<double> eigen(N);
	
        Vector data(N);
	for(size_t i = 0; i<N; ++i){
            double xpos =L*0.5-double(i)*L/double(N);
		data[i]=1./sqrt(2*pi)*exp(-xpos*xpos/2.);
	}
	matrix<double> H(N,N);
	for(size_t i=0; i<H.size1();++i){
            double q = pi*40/L;
		double xpos =L*0.5-double(i)*L/double(N);
		double sinx =sin(q*xpos);
		for(size_t j=0; j<H.size2();++j){
		
H(i,j)=(i==j?1/(deltax*deltax)+g*data[i]*data[i]+xpos*xpos*0.5+q*q*0.5*sinx*sinx
				:((i==j-1) || (i==j+1)? -0.5/(deltax*deltax):0) );
		}
	}
	//std::cout<<H;
	
	char jobz = 'V';
	char uplo = 'U';
	int n = H.size1();
	//double* A = new double[];
	//A=m;
	int lda = n;
	int lwork = 3*n-1;
	double* work = new double[lwork];	
	int info=2313;
	
	
	dsyev_(&jobz, &uplo, &n, &H(0,0), &lda, &eigen[0], work, &lwork, &info);
	
	//copy(eigen.begin(), eigen.end(), ostream_iterator<double>(cout,"\n"));
	//cout<<"&"<<endl;	
	for(int j=0; j<3;++j){
		for(int i=0; i<H.size1();++i)
			cout<<H(j,i)<<endl;
		cout<<"&"<<endl;
		
	}
		
	return 0;
}
