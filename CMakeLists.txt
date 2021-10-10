cmake_minimum_required(VERSION 3.20)
project(Pliki_ródowe_20211009)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(
        lab_no.h
        main.cpp
        matrix.cpp
        matrix.h
        ode_solver.cpp
        ode_solver.h
        opt_alg.cpp
        opt_alg.h
        solution.cpp
        solution.h)
