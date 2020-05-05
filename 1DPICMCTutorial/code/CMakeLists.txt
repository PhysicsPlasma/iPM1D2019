list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake-modules")
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})

file(GLOB_RECURSE SRCS base/* mc/* pic/* diagnostics/* solver/*.c *.F90)
add_executable(${PROJECT_NAME} ${SRCS})

include(FetchContent)

FetchContent_Declare(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG master  # whatever desired version is
)

FetchContent_MakeAvailable(h5fortran_proj)

# ------------------------------------------------------
# whatever your program is
target_link_libraries(${PROJECT_NAME} h5fortran::h5fortran)
target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/1DPICMCTutorial)
