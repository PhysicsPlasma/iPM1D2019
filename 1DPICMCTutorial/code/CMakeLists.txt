# 设置项目所需cmake的最低版本
cmake_minimum_required(VERSION 3.0.0)
# 项目名称，也是最后生成的可执行文件名
project(1DPICMIC)
# 启用对Fortran语言的支持
enable_language(Fortran)
# find_package(MPI)
# 设定编译组态，RelWithDebInfo指优化运行速度同时产生debug信息
set(CMAKE_BUILD_TYPE RelWithDebInfo)
#将当前项目里指定目录下的所有文件(包含子目录下的文件)添加到SRCS这个变量里
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} base/* diagnostics/* mc/* pic/* *.F90)
# 将SRCS编译为mylib这个动态链接库
add_library(mylib SHARED ${SRCS})
# 使用${SRCS}中的源文件生成可执行文件
add_executable(1DPICMIC ${SRCS})
# 让myprogram依赖mylib
add_dependencies(1DPICMIC mylib)
# 链接mylib
target_link_libraries(1DPICMIC mylib)
include(FetchContent)

FetchContent_Declare(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG master  # whatever desired version is
)

FetchContent_MakeAvailable(h5fortran_proj)

# ------------------------------------------------------
# whatever your program is
target_link_libraries(1DPICMIC h5fortran::h5fortran)
