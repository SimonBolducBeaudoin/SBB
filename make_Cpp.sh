#!/bin/bash

#make -C Cpp/Autocorrelation/ -j4
#make -C Cpp/Benchmark_test/ -j4
#make -C Cpp/Bispectrum/ -j4
#make -C Cpp/Convolution/ -j4
#make -C Cpp/FFT/ -j4
make -C Cpp/FFTW_extra/ -j4 
make -C Cpp/Math_extra/ -j4 
make -C Cpp/Omp_extra/ -j4 
make -C Cpp/Moments_cumulants/ -j4 
make -C Cpp/Numerical_integration/ -j4 
make -C Cpp/Special_functions/ -j4 
make -C Cpp/Windowing/ -j4 
make -C Cpp/aCorrs-OTF/ -j4 
make -C Cpp/Time_quadratures/ -j4 
make -C Cpp/Histograms/ -j4 




