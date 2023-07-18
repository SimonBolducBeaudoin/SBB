#!/bin/bash

#make -C AFAP/Autocorrelation/ -j4
#make -C AFAP/Benchmark_test/ -j4
#make -C AFAP/Bispectrum/ -j4
#make -C AFAP/Convolution/ -j4
#make -C AFAP/FFT/ -j4
make -C AFAP/FFTW_extra/ -j4 
make -C AFAP/Math_extra/ -j4 
make -C AFAP/Omp_extra/ -j4 
make -C AFAP/Moments_cumulants/ -j4 
make -C AFAP/Numerical_integration/ -j4 
make -C AFAP/Special_functions/ -j4 
make -C AFAP/Windowing/ -j4 
make -C AFAP/aCorrs-OTF/ -j4 
make -C AFAP/Time_quadratures/ -j4 
make -C AFAP/Histograms/ -j4 




