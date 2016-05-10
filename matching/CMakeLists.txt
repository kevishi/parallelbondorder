set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -g")
cmake_minimum_required (VERSION 2.6)
project (pmatch)

install(FILES
  "${PROJECT_BINARY_DIR}/lapack-3.6.0/lapack_build.cmake"
 DESTINATION "./" COMPONENT dev)
include_directories(lapack-3.6.0/LAPACKE/include lapack-3.6.0/LAPACKE/example lapack-3.6.0/LAPACKE/utils lapack-3.6.0/LAPACKE/src eigen-eigen-07105f7124f9/Eigen)
file (GLOB source src/*.cpp)
add_executable(pmatch ${source})
