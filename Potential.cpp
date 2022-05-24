#include <iostream>
#include <vector> 
#include <algorithm>
#include <functional>
#include <cmath>
#include <climits>
#include <iterator>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

int main (){
    const double pi = boost::math::constants::pi<double>();
    double L = 10.;
    int N = 1024;
    double deltax = L/double(N);
    double g=1;
    
    vector<double> pot(N);
    for(size_t i=0; i<pot.size();++i){
        double q = pi*40/L;
        double xpos =L*0.5-double(i)*L/double(N);
        double sinx =sin(q*xpos);
        pot[i]=xpos*xpos*0.5+q*q/2.*sinx*sinx;
    }
    
    
    copy(pot.begin(), pot.end(), ostream_iterator<double>(cout, "\n") );
    
    return 0;
    
}
                
                
        