#ifndef LEGENDRE_H
#define LEGENDRE_H

template<int N>
struct LegendrePolynomial{
    inline static double eval(double x){
        return ((2.*(N-1)+1)*x*LegendrePolynomial<N-1>::eval(x)-(N-1)*LegendrePolynomial<N-2>::eval(x))/(N);
        }
};

struct LegendrePolynomial<1>{
    inline static double eval(double x){
        return x;
    }
};

struct LegendrePolynomial<0>{
    inline static double eval(double x){
        return 1.;
    }
};

inline double LP(double x, int N){
    return ((2.*(N-1)+1)*x*((N-1==1)?x:LP(x, N-1)))-(N-1)*((N-2==0)?1.:LP(x, N-2))/(N);
}

    
#endif