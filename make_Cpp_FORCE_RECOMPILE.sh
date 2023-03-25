#!/bin/bash

#make -C Cpp/Autocorrelation/ -j4 -B
#make -C Cpp/Benchmark_test/ -j4 -B
#make -C Cpp/Bispectrum/ -j4 -B
#make -C Cpp/Convolution/ -j4 -B
#make -C Cpp/FFT/ -j4 -B
make -C Cpp/FFTW_extra/ -j4 -B
make -C Cpp/Math_extra/ -j4 -B
make -C Cpp/Omp_extra/ -j4 -B
make -C Cpp/Moments_cumulants/ -j4 -B
make -C Cpp/Numerical_integration/ -j4 -B
make -C Cpp/Scoped_timer/ -j4 -B
make -C Cpp/Special_functions/ -j4 -B
make -C Cpp/Windowing/ -j4 -B
make -C Cpp/aCorrs-OTF/ -j4 -B
make -C Cpp/Time_quadratures/ -j4 -B
make -C Cpp/Histograms/ -j4 -B




