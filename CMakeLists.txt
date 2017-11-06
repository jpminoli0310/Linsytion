cmake_minimum_required(VERSION 3.0)
project(Linsition Fortran)

set(example Example1)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

SET(EXECUTABLE_OUTPUT_PATH ../run/)
add_executable(${example}.out ../Examples/${example}.f90 ../src/Linsytion.f90)
target_link_libraries(${example}.out ${OpenMP_LIBS})
