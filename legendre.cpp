#include<iostream>
#include"headers/legendre.h"
#include<vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

extern "C" void dstev_(char* jobz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info); 
typedef boost::numeric::ublas::matrix<double> Matrix;

int main(){
    char jobz='N';
    int N = 10;
    const int legdeg = 10;
    Matrix Rmat(N, N);
    int ldz = N;
    int info;
    double *work = new double[2*N-2];
    vector<double> d(N, 0);
    vector<double> e(N);
    for(int i = 0; i<N; ++i)
        e[i]=((i+2)-1.)/sqrt((2.*(i+2.)-3.)*(2.*(i+2)-1));
    dstev_(&jobz, &N, &d[0], &e[0], &Rmat(0,0), &ldz, work, &info);
    for(int i = 0; i<N; ++i)
        cout<<d[i]<<'\t'<<LegendrePolynomial<legdeg>::eval(d[i])<<endl;
}