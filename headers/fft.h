#ifndef FFT_H
#define FFT_H

#include <vector>
#include <iterator>
#include <complex>
#include "fftw3.h"


namespace fft {

  struct Forward  {
    static const int d = FFTW_FORWARD;
  };

  struct Backward {
    static const int d = FFTW_BACKWARD;
  };

  template<class Direction>
  struct TransformDirection{
    static const int direction = Direction::d;
  };


  typedef TransformDirection<Forward>  forward;
  typedef TransformDirection<Backward> backward;

  struct Measure {
    static const unsigned flag = FFTW_MEASURE;
  };

  struct Estimate {
    static const unsigned flag = FFTW_ESTIMATE;
  };

  template<
    class DirectionPolicy,
    class PlanPolicy,
    class VECTOR
    >
  class FTImpl{
    public:
      FTImpl(const VECTOR& input, int dimension);
      //copy constructor does suit so far...
      ~FTImpl();

      std::vector<std::complex<double> > exec();



    private:
      int N_, dim_;
      fftw_complex* in;
      fftw_complex* out;
      std::vector<std::complex<double> > cop_;
      std::vector<std::complex<double> >::iterator it;
      fftw_plan p;
};

///////////////////////////// Implementation starts here... ///////////////////////////////////////////////////

template<class DirectionPolicy, class PlanPolicy, class VECTOR>
FTImpl<DirectionPolicy, PlanPolicy, VECTOR>::FTImpl(const VECTOR& input, int dimension):	
N_(input.size()), cop_(N_) , dim_(dimension){
	//std::copy(input.begin(), input.end(), cop_.begin());
    in=static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N_) );
    out=static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N_) );
    for(int i = 0; i<N_; ++i){
        in[i][0]=real(input[i]);
        in[i][1]=imag(input[i]);
    }
        
        p=fftw_plan_dft(dimension, &N_, in, out, DirectionPolicy::direction, PlanPolicy::flag);

}

template<class DirectionPolicy, class PlanPolicy, class VECTOR>
FTImpl<DirectionPolicy, PlanPolicy, VECTOR>::~FTImpl(){
    fftw_free(in);
    fftw_free(out);
	fftw_destroy_plan(p);
	fftw_cleanup();
}

template<class DirectionPolicy, class PlanPolicy, class VECTOR>
std::vector<std::complex<double> > FTImpl<DirectionPolicy, PlanPolicy, VECTOR>::exec(){
	fftw_execute(p);
	for(int i = 0; i<N_; ++i)
		cop_[i] = std::complex<double>(out[i][0], out[i][1]);						
	return cop_;		
}

}//namespace fft

template<class DirectionPolicy, class PlanPolicy, class VECTOR>
inline std::vector<std::complex<double> > FourierTransform(VECTOR& input, int dim = 1){
	fft::FTImpl<DirectionPolicy, PlanPolicy, VECTOR> ft(input, dim); 
	std::vector<std::complex<double> > ret =  ft.exec();
	return ret;
}


#endif
