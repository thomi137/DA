#include<boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include<vector>
#include<iostream>
#include<algorithm>

//using namespace std;
using namespace boost::numeric::ublas;

int main(){
	
	vector<double> v(3);
	vector<double> v1(3);
		for(int i = 0; i<v.size(); ++i){
			v(i)=i;
			v1(i)=i;
		}
		
	std::cout<<outer_prod(v,v1);	
	vector_range<vector<double> > vr(v, range(1,2) );
	vr(1)=1;
	std::cout<<vr;
	vector<double>::iterator it = v.begin();
	std::vector<double> vec(3); 
	//v=prod(m, v);
	std::copy(v.begin(), v.end(), vec.begin());	
	//std::cout<<3*vr;				

	std::copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout,"\n") );
	
	return 0;
}			
		
