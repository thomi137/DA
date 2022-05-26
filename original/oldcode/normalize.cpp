#include<vector>
#include<complex>
#include<iostream>

using namespace std;

template <class InputIterator1, class T>
	inline T norm_of_vec( InputIterator1 first, InputIterator1 last, T init){
	for(;first != last;++first)
		init+=norm(*first);
	return init;
	}
			

int main(){
	typedef complex<double> ct;
	int N = 1024;
	vector<ct> vec(N, 1./sqrt(2*M_PI) );
	
	cout<<2*M_PI/N*norm_of_vec(vec.begin(), vec.end(), 0.)<<endl;
	
	return 0;
	
}
