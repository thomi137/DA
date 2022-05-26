#include<iostream>
#include<vector>
#include<iterator>
#include<functional>
#include<algorithm>
#include"fftw3.h"
#include<complex>

using namespace std;

int main(){
	typedef complex<double> cmpl;
	cmpl I(0,1);
	vector<cmpl> vec1(14, I);
	vector<cmpl> vec2(14, 5.);

	copy(vec1.begin(), vec1.end(), ostream_iterator<cmpl>(cout,"\n")
			);
	
	transform(vec1.begin(), vec1.end(), vec2.begin(), vec1.begin(),
			multiplies<cmpl>());
	
	copy(vec1.begin(), vec1.end(), ostream_iterator<cmpl>(cout,"\n")
			);

	transform(vec2.begin(), vec2.end(), vec2.begin(),
			bind2nd(divides<cmpl>(),5)
			);
	
	copy(vec2.begin(), vec2.end(), ostream_iterator<cmpl>(cout,"\n")
			);
	
	}
	
