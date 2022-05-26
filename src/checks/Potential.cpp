#include <iostream>
#include <vector> 
#include <algorithm>
#include <cmath>
#include <iterator>

#include <boost/math/constants/constants.hpp>

// compile with, e.g.
// c++ --std c++11 -I /usr/local/include/boost Potential.cpp -o Potential
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
        cout << xpos << ',' << pot[i] << endl;
    }
    return 0;
}
                
                
        