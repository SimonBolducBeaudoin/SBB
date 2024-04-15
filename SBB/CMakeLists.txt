cmake_minimum_required(VERSION 3.15)
project(AutoCorr) 

if(MSVC)
    # ONLY CUDA CODE ON MSVC NOTHING ELSE
    add_subdirectory(FFT)
elseif(CYGWIN)
    # NO CUDA AVAILABLE ON CYGWIN UNFORTUNATELY
    add_subdirectory(FFT)
    #add_subdirectory(Convolution) #tbd
    add_subdirectory(Multi_array)
    add_subdirectory(Math_extra)
    add_subdirectory(Omp_extra)
    add_subdirectory(AutoCorr)
    add_subdirectory(Histograms)
    add_subdirectory(Time_quadratures)
else() #unix
    add_subdirectory(FFT)
    #add_subdirectory(Convolution) #tbd
    add_subdirectory(Multi_array)
    add_subdirectory(Math_extra)
    add_subdirectory(Omp_extra)
    add_subdirectory(AutoCorr)
    add_subdirectory(Histograms)
    add_subdirectory(Time_quadratures)
endif()

