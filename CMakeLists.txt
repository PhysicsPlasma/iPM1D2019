cmake_minimum_required(VERSION 3.10.0)
project(1DPICMC LANGUAGES C Fortran)
add_subdirectory(1DPICMCTutorial)
file(COPY
        ${PROJECT_SOURCE_DIR}/1DPICMCTutorial/gas
        ${PROJECT_SOURCE_DIR}/1DPICMCTutorial/input
    DESTINATION
        ${PROJECT_SOURCE_DIR}/run)