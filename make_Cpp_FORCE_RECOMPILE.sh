#!/bin/bash

#make -C AFAP/Autocorrelation/ -j4 -B
#make -C AFAP/Benchmark_test/ -j4 -B
#make -C AFAP/Bispectrum/ -j4 -B
#make -C AFAP/Convolution/ -j4 -B
#make -C AFAP/FFT/ -j4 -B
make -C AFAP/FFTW_extra/ -j4 -B
make -C AFAP/Math_extra/ -j4 -B
make -C AFAP/Omp_extra/ -j4 -B
make -C AFAP/Moments_cumulants/ -j4 -B
make -C AFAP/Numerical_integration/ -j4 -B
make -C AFAP/Scoped_timer/ -j4 -B
make -C AFAP/Special_functions/ -j4 -B
make -C AFAP/Windowing/ -j4 -B
make -C AFAP/aCorrs-OTF/ -j4 -B
make -C AFAP/Time_quadratures/ -j4 -B
make -C AFAP/Histograms/ -j4 -B




