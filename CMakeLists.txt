cmake_minimum_required(VERSION 3.20)
project(optimization)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(optimization
        src/lab_no.h
        src/main.cpp
        src/matrix.cpp
        src/matrix.h
        src/ode_solver.cpp
        src/ode_solver.h
        src/opt_alg.cpp
        src/opt_alg.h
        src/solution.cpp
        src/solution.h)
