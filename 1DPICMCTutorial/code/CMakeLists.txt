file(GLOB_RECURSE SRCS base/* diagnostics/* mc/* pic/* *.F90 *.f90)
add_executable(1DPICMCStartup ${SRCS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/run)