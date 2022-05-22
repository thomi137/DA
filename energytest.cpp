#include "./headers/bec.h"

using namespace std;

int main(){
		
	//testing the energy function of bec.h
	int N=64;
	vector<complex<double> > gs(N, 1./sqrt(2.*M_PI) );
	
	cout<<bec::energy(gs, 0.);
	
	return 0;
	
}
