
# SBB
A python module for all the python code I'm devellopping for my lab experiements.

# Installation
## Before installing
- Make sure Python and pip are installed on your system
- On windows the use of cygwin is recommended

# Dependencies
- Other repositories from me
  
They are place in ../Cpp a repository parallel to SBB/

# Dependencies
    All homebrewed libraries are imported using global imports "#include <library.h> ".
    This means that the library must eather be installed in your environnment's path or that it must be included during compilation and linking. 
    You can edit the "CmakeList.txt" to properly include homebrewed libraries.
    - Homebrewed libraries (available on my github : https://github.com/SimonBolducBeaudoin)
        - aCorrs-OTF
        - Convolution
        - FFT
        - FFTW_extra
        - Histograms
        - Math_extra
        - Moments_cumulants
        - Numerical_integration
        - Omp_extra
        - Scoped_timer
        - Special_functions
        - Time_quadratures
        - Windowing
    - Other dependencies<
        - openmp
        - mprf (mpreal)
        - fftw3
            - Can be installed using your package manager.
        - pybind11
    Pybind11 can be installed using you're python package manager (conda(anaconda env),pip,pacman,...).
    
# Building, compiling and Installing
    - Edit config.cmake for your machine (If you are compiling in a different envionnment than your python installation) so that pybind11 can be detected and used.
    - Unix environnment
		- cmake -S . -B ./build && cmake --build build/ && cmake --install build/
    - Crosscompiling to for windows (Cygwin)
		Pass the toolchain to cmake to build the project (the rest is the same).
        - cmake -S . -B ./build -DCMAKE_TOOLCHAIN_FILE=./CMakeConfigs/mingw_toolchain.cmake && cmake --build build/ && cmake --install build/
    
	
	
# Removing the build directory
    cmake doesn't offer a built-in solution. 
    Best solution is to use rm.
    - rm -R -f build/
    
#  Cleaning the compilation's output (equivalent to make clean)
    - cmake --build build/ --target clean
 
