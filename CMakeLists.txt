cmake_minimum_required(VERSION 3.15)
project(SBB) 

message(STATUS "Building AFAP")
add_subdirectory(AFAP)

message(STATUS "Installing AFAP")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/AFAP/aCorrs-OTF/obj/
		DESTINATION ${PROJECT_SOURCE_DIR}/SBB/AutoCorr/
		FILES_MATCHING PATTERN "*.pyd;*.py")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/AFAP/Time_quadratures/obj/
		DESTINATION ${PROJECT_SOURCE_DIR}/SBB/Time_quadratures/
		FILES_MATCHING PATTERN "*.pyd")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/AFAP/Histograms/obj/
		DESTINATION ${PROJECT_SOURCE_DIR}/SBB/Histograms/
		FILES_MATCHING PATTERN "*.pyd")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/AFAP/Moments_cumulants/obj/
		DESTINATION ${PROJECT_SOURCE_DIR}/SBB/Moments_cumulants/
		FILES_MATCHING PATTERN "*.pyd")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/AFAP/Omp_extra/obj/
		DESTINATION ${PROJECT_SOURCE_DIR}/SBB/Omp_extra/
		FILES_MATCHING PATTERN "*.pyd")


        
