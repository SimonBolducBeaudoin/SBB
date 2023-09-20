cmake_minimum_required(VERSION 3.15)
project(AutoCorr) 

message(STATUS "Building SBB")
add_subdirectory(Math_extra)
add_subdirectory(Omp_extra)
add_subdirectory(AutoCorr)
add_subdirectory(Histograms)


