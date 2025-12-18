cmake_minimum_required(VERSION 3.31)
project(LabNumericalMethods4)

set(CMAKE_CXX_STANDARD 20)

add_executable(LabNumericalMethods4 Task1.cpp

        Methods/Eiler.h

        Methods/RungeKutta4.h)
