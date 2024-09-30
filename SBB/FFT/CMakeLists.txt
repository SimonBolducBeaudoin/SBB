# Should be compatible with unix and 
 cmake_minimum_required(VERSION 3.15)
project(FFT) # Needs to be ontop
if(MSVC)
    # ONLY CUDA CODE ON MSVC NOTHING ELSE
    add_subdirectory(FFT_CUDA)
elseif(CYGWIN)
    # NO CUDA AVAILABLE ON CYGWIN UNFORTUNATELY
    add_subdirectory(FFTW_extra)
    add_subdirectory(DFT)
else() #unix
    add_subdirectory(FFT_CUDA)
    add_subdirectory(FFTW_extra)
    add_subdirectory(DFT)
endif()
