
#include"bec.h"

using namespace std;


int main(){
  int N = 128;
  
  complex<double> I = bec::I, tempo;
  vector<complex<double> > vec(N);
  for(int i = 0; i < N; ++i){
    vec[i]=complex<double>(sin(2.*M_PI/double(N)*i));
  }
  for(int i =0; i<vec.size();++i)
    cout<<real(vec[i])<<endl;
  bec::GP_Hamiltonian(vec, 0.1);
  cout<<"&"<<endl;
  for(int i =0; i<vec.size();++i)
     cout<<real(vec[i])<<endl;
 



 //  complex<double> e=bec::energy(vec, 0.1);
			       
//   cout<<e<<endl;
  
}
