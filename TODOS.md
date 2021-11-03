## Conceptuals
- Problems
    + (16-12-19) With current calculation methods the time of claculation 
			(time_quads + hist) is of the same order (but still higher ) than 
			acquisition. Its only going to get worst when adding multiple 
			numecrical filters.
- Questions
	+ (16-12-19) What cost the most in computation time in a 
			aquasition + time_quad (single of multi-kernels) + histograms ? 
			And what could potentially be reduced ?
- Ideas ?
    + (16-12-19) Using SIMD operations to speedup multi-kernels operations 
			on the same data vector ?
    + (16-12-19) Using CUDA for even faster FFTs (FFTWs) ?
## Paperworks
- Report
    - (12-12-19) Demonstration on analytical cases
    - (16-12-19) Revive my latex thesis
        - (10-01-20) Done
    - Start theory for linear amplifier in time
        - (10-01-20) Done
    - (15-01-20) Multivariate culumant for the blackbody experiment
        - (16-01-20) opened
        - () writing theory on multivariate cumulants 
    - (16-12-19) Summery of the project ?
- Oral presentation
    - (23-01-20) Group meeting 14 fev
        + Start with narrow bandwidth relations for quadratures
        + Then Reulet transforms
            + Define second quantization basis
            + Example in frequency representation
            + Example in time and wavelets representation
            + Link between \hat[a] , \hat[v], \hat[X]_{Theta} continuum rep.
        + Demo for black body statistics in wavelets rep.
        + How to get this statistic experimentally
            + Microwave experiment
                + Cumulant properties
                + Monovariate cumulant
                    + up to 4th order
                + Bivariate cumulant
                + Application on the black body experiment
                    + Setup
            + Optic experiment
                + Measurment directly in Fock basis
                + Narrow bandwidth
                + Setup ?
        + Measure of state representable in 2nd quantification basis
            + Large bandwidth squeezed state
            + Finding the maximaly squeezed basis
            + Proposition of approach ?
        + Other state ?

## Workflow
- Git stuff
    - (16-12-19) Find how to properly update my project without always having the "detached-HEAD" bug.
	 + Stop using git submodules 
    - Automate git commits (and push ?) before making actual measurments ?
    - Git doesn't distinguish between Caption.txt and caption.txt (it anoying)
- Python stuff
    - Redesign the Python directory ?
    - Learn pycharm ?
    - Learn spyder ?
- Notebook 
    - Use jupiter to summurized experiements (after or during experiement ?)
        + dont use jupyter but use markdown inside code ? and combine with markdown viewer inside notepad++
- Compiling stuff
    - Find a more intuitive, more flexible solution to replace gnu-make
- C++ 
    - Explicitly write down my directory/file conventions
- gvim 
    - look it up and test
- premake
    - A language to manage the generation of makefile
- markdown
    - Is there a for of mardown that would help to organize c++ text file and python text file ?

## Python

### Proof of concept

### Validations
- moment_k_float centered_moment_k_float and cumulant_k in Python
    + (13-01-20) OPENED
    + For Histogram1D
    + Validated up to 3rd order using numpy.random.normal
    + (13-01-20) CLOSED

### Benchmarks
- Re-make a histogram benchmark and save it
    * (13-12-19) OPENED
    + (16-12-19) 1D done and saved under Benchmark/Results;
    + (16-12-19) 2D done and saved under Benchmark/Results;
    * (16-12-19) CLOSED
- (15-12-19) Benchmark the measure + quad + hist process

### Routines
- Guzik V measurment 
    + (12-12-19) DONE
- Calib Gain routine
    + (12-12-19) Manual version DONE
- Quadrature measurement
    + (12-12-19) DONE
- Quadrature Histograms
    * (12-12-19) OPENED
    + 1D
    + 2D
    * (12-12-19) CLOSED 
- Quadrature Histograms with multiple measures
    * (12-12-19) OPENED
    * (13-12-19) CLOSED
- Parrallel aquisition//Traitement
    * (13-12-19) OPENED
    + using python threads ( Théo. ~1.6 fois plus rapide)
- Temperature control
    * (16-12-19) OPENED
    + (13-12-19) Formation avec Ed
    + Using class inheritance from pyHegel's existing Lakeshore370
    + (20-12-19) Test Virtual instrument "temperature_control"
    + Proper PyHegel inherintance DONE
    + Round of test to calibrate tables
        + (20-12-19) Round 1 Done
        + (10-01-20) Round 2
            ++ Tested for 7 to 200 mK. 
            ++ Tables : temperature_range, heater_range , PID, etc... are tested but could be calibrated furthermore
    + Auto stabilization (10-01-20)       
    * CLOSED
    * BUG
    + At object creation :""
    + At object destruction :""
    * TODOS
    + Change class name to "lakeshore_370_logical"
    + Change method name to table_set_and_get
    + Add monitor option for live assesment
    + Calibrate up to 800 mk

### Experiments

## C/C++
- Omp_extra
	- (28-01-20) Add a function that can be called to check if the user asks for more threads
		then are possibly available and gives the max amount of threads when it aksed to much.
- Math_extra
    - (13-01-20)Created
- Moment_cumulants
    - TASK_01 Adding bivariate case to SM-Moments_cumulants
- Histogram class
    - (12-12-19) properly accomodate int 8 and int16 in 1D and 2D
    - (12-12-19) abscisse function for python
        +  Done for double type 
    - (12-12-19) Add a memory_size function
    - (13-12-19) BUG : 1D histograms.sum() is not equal to the size of the input
                double as inputs
        + (16-12-19) Fixed : #pragma reduction makes copies of the reduction variable for each thread 
                            but when this variable is also a member attribute of the class and is called inside a function 
                            called by the child thread the member attribute will be called instead of the copy made for the thread. 
                            This can lead to race conditions in the code.
        + (16-01-20) This bug was supposed to be fixed but it wasn't now it should be fixed.
    - (13-01-20) Pybind11 wrapper for 1D moment_k and cumulant
    - TASK_02 Adding first order moment and cumulant to 2D histograms (OPENED)
	- (28-01-20) Return an error message when user tries to make 2D hitogram with nofbins bigger than 2**12
	- Optimize to minimize for cache misses
		
- time_quadratures class
    - (12-12-19) Transfert Fresnel Sine and Cos to a special function submodule
    - (12-12-19) Add a memory_size function 
    - (16-12-19) Convert time_quadrature class into python callable function and
            create SM-Blackbody which will contain a combined code 
            specialized for the blackbody experiment.
	- (28-01-20) Add custom filter inputs 
- Special function submodule
    - (16-12-19) Fonction de filtrage réelle et gaussienne (autre forme ?) autour 
            deux ou trois fréquences.
- Blackbody submodule
    - (16-12-19) Create SM-Blackbody from time_quadrature
    - (16-12-19) Accumuler l'histogram a mesure que les quadratures
            sont calculées.
    - (16-12-19) Pré-calculer les kernels pour plusieurs filtres numériques.
    - (16-12-19) Minimiser le coût associer à l'utilisation de plusieur kernels et filtres.
- Multi_array submodule
	-
- Time_domain library
    - (12-12-19) Make the makefile

# Lab book

## 24-01-20

TASK_01 Adding bivariate case to SM-Moments_cumulants
- General bivariate moments implemented DONE(28-01-20)
- Only cumulant  c_ii c_jj c_iijj implemented DONE(28-01-20)

TASK_02 Adding first order moment and cumulant to 2D histograms
- Pybind11 wrapper working (see TASK_01) DONE(28-01-20)
- Validation/Moments_cumulants_2.py creatre
	- monovariate moments of 2D histograms were validated with scipy.stat.moment

# 10-02-20
- Only l_data, l_valid, l_full and bin depth of 1D histograms should be uint64_t other should be uint
- 

# 24-02-20
- Removed git submodule
- Added a feature to Histograms
	- Cap the number of thread depending on machine
- Make it so that TimeQuad object return only valid part to Python
	- Done and tested
- Prepare routine PQvsVDC_0
	- Aquiering  p and q using multiple filters
	- Accumulating the results in 2D Histogram leeping track of experimental parameters

# 26-02-20
- Re-read TDQO-V1
- Print plan d'etude
- Prepare 2-Slides on TDQO measurments
- Open a long beamer presentation

# 29-02-20
- Determine  how many ADC there are using histograms
    - 2 ADCs with different response

 
# 02-03-20
- Bug fix in TimeQuad
    - (Fixed) Loopping on TimeQuad.execute(data) was always accumulating in ps and qs (they kept on growing)
        - Fixed by reseting ps and qs to 0 in between execution (only overlapping part between threads)
          and using = instead of += when possible.
    - (Fixed) Was not behaving as indtended when 1 filter was given...
- Wrote Validation script TimeQuad_7 and 8 and PQsVsVdc_0
- Bug fix in Multi_array
    - (Temporally fixed) Multi_array's index overflowed when array is big.
- Launch an overnight measurment of expval n and n^2 vs Vdc
   - k
- Write a quick implementation for bistrpectrum

# Futur
- Write a Guzik logical device for pyHegel
- Same for yoko
- Update Timing labrary to use high precision time (see: fftw  documentation on timing)
- Time_quad can probably by optimized 
    - SIMD instruction could be used
    - Using fftw's multi-dimension fft to manage fft's in bulk ?
    - Optimization for FFT size and l_data not trivial
- Checks if n_threads checks works properly for histogram class
- Prepare routine PQvsVDC_0
	- Calculating cumulants from those
	- Display expval of n and n^2 from thos data
- Write down the clean demonstration for 2 coherents sources
- Histograms: BUG FIX to implement
    - histograms dont properly manage inputs of the for numpy_array[s,S,Step] because it only uses ptr and size
        - Fix : Implement Multi_array<DataType,1> data 
            that is initialized empty and 
            that then uses the asignment operator to inherit the properties of a np array
            - Multi_array<T,1> operator= is written (no test done)
            - Next step : add member Multi_array<DataType,1> data in histogram
            - Test : it should avoir unnecessary copies
            - Implement for Dim = 2 and 3.
    - Also Replace all memory manegement in hitograms for Multi_arrays
- Bug in Multi_array : there can be some overfloa in some condition. 
    - Need to add check at creation
    - Need to optimize for speed

# OPENNED BUGS 
- Pybind11 doesn't work with uint64_t 
	- This coul potentially cause problems with l_data, l_valid and l_full

