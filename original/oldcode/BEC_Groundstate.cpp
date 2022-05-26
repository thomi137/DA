#include <iostream>
#include <vector> 
#include <algorithm>
#include <functional>
#include <cmath>
#include <climits>
#include <iterator>
#include <ietl/lanczos.h>
#include <ietl/vectorspace.h>
#include <ietl/interface/ublas.h>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric::ublas;


typedef complex<double> cmpl;
typedef boost::numeric::ublas::vector<double> Vector;
typedef std::vector<double> StlVec;
typedef boost::lagged_fibonacci607 Gen;
typedef ietl::vectorspace<Vector> Vecspace;
  	    
template<class VecType>
class Hamiltonian{
	public:
		typedef VecType Vector;	
		typedef typename VecType::value_type type;
		typedef matrix<type> MatType;	
		Hamiltonian(int N, double L, double g, VecType psi):
			N_(N), L_(L), g_(g), psi_(psi), V_(N), Lap_(N,N){
		deltax = double(N)/L;
		for(unsigned int i = 0; i<V_.size(); ++i){
			double xpos =L*0.5-double(i)*deltax;
			V_(i)=g_*psi_(i)*psi_(i)+xpos*xpos/2.;
			}
			
		for(unsigned int i = 0; i<Lap_.size1();++i){
			for (unsigned int j=0; j<Lap_.size2();++j){
				Lap_(i,j) = (i==j? -1/(deltax*deltax)+V_(i):
			                    ((i==j-1) || (i==j+1)? 0.5/(deltax*deltax):0.) );
				}           				    					
			}
		}	
	 
	void mult (const Vector& x, Vector& y) const {
			y=prod(Lap_, x);
			y(0)=y(N_-1)=0.;
		}				
	private:
	Vector V_, psi_;
	MatType Lap_;
	int N_;
	double L_, g_, deltax;
};


namespace ietl{	
	inline void mult(const Hamiltonian<Vector>& h, const Vector& x, Vector& y){
		h.mult(x,y);
	}
}//namespace itl



int main (){
	const double pi = 4*atan(1.0);
	double L = 50.;
	int N = 32;
	double deltax = L/double(N);
	double g=0.;
	int n_lowest_eigval=2;
	int max_iter = 50000;
	
	Gen mygen;
	
	double rel_tol = 500.*(numeric_limits<double>::epsilon());
	double abs_tol = pow(numeric_limits<double>::epsilon(), 2./3.);
	
	std::vector<double> eigen;
	std::vector<Vector> eigenv;
	std::vector<double> err;
	std::vector<int> multiplicity;
	
        Vector data(N);
	for(size_t i = 0; i<N; ++i){
            double xpos =L*0.5-double(i)*L/double(N);
		data[i]=1./sqrt(2*pi)*exp(-xpos*xpos/2.);
	}
	
	copy(data.begin(), data.end(), ostream_iterator<double>(cout,"\n"));
	cout<<"&"<<endl;	
		
	Vecspace vec(N);
	Hamiltonian<Vector> ham(N, L, g, data);
	ietl::lanczos<Hamiltonian<Vector>, Vecspace> lanczos(ham, vec);
	ietl::lanczos_iteration_nlowest<double> iter(max_iter, n_lowest_eigval, rel_tol, abs_tol);

	try{
	   lanczos.calculate_eigenvalues(iter, mygen);
	   eigen=lanczos.eigenvalues();
	   err=lanczos.errors();
	   multiplicity = lanczos.multiplicities();
  	 }
	 catch(std::runtime_error& e){
		 cout<<e.what()<<endl;
	 }
	
	 
	 ietl::Info<double> info;
	 
	 try{
	    lanczos.eigenvectors(start, end, back_inserter(eigenv), info, mygen);
    	}
    	catch(std::runtime_error& e){
		cout<<e.what()<<endl;
	}
	
	for(std::vector<Vector>::iterator it = eigenv.begin(); it != eigenv.end(); ++it){
		copy(it->begin(), it->end(), ostream_iterator<double>(cout, "\n") );
		cout<<"&"<<endl;
	}
	
	return 0;
}
