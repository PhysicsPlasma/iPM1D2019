cmake_minimum_required(VERSION 3.12.0)
set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)
project(iPM1D2019  LANGUAGES C Fortran)
set(CMAKE_BUILD_TYPE "Release")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coarray=shared -coarray-num-images=4")
add_compile_options(-heap-arrays)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,stack-size=1000000")

file(REMOVE_RECURSE
        ${PROJECT_SOURCE_DIR}/run
        ${PROJECT_SOURCE_DIR}/build
        ${PROJECT_SOURCE_DIR}/run/restart)

file(GLOB_RECURSE SRCS "1DPICMCTutorial/code/base/*"
                    "1DPICMCTutorial/code/diagnostics/*"
                    "1DPICMCTutorial/code/mc/*"
                    "1DPICMCTutorial/code/pic/*"
                    "1DPICMCTutorial/code/*.F90")

add_executable(${PROJECT_NAME} ${SRCS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/run)

file(COPY
         ${PROJECT_SOURCE_DIR}/1DPICMCTutorial/gas
         ${PROJECT_SOURCE_DIR}/1DPICMCTutorial/input
  DESTINATION
        ${PROJECT_SOURCE_DIR}/run)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/run/restart)
