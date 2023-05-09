#!/bin/bash

# This bash file is used to build everything for the SBB module

#make -C Cpp/Autocorrelation/ clean
#make -C Cpp/Bispectrum/ clean
#make -C Cpp/Convolution/ clean
#make -C Cpp/FFT/ clean
make -C Cpp/FFTW_extra/ clean 
make -C Cpp/Math_extra/ clean 
make -C Cpp/Omp_extra/ clean 
make -C Cpp/Moments_cumulants/ clean 
make -C Cpp/MPFR_extra/ clean 
make -C Cpp/Multi_array/ clean 
make -C Cpp/Numerical_integration/ clean 
make -C Cpp/Special_functions/ clean 
make -C Cpp/Windowing/ clean 
make -C Cpp/aCorrs-OTF/ clean 
make -C Cpp/Time_quadratures/ clean 
make -C Cpp/Histograms/ clean 





