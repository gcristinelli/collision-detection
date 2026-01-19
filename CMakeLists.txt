cmake_minimum_required(VERSION 4.0)
project(Assignment)

set(CMAKE_CXX_STANDARD 20)

add_executable(Assignment main.cpp
        Motion.cpp
        Motion.h
        Output.cpp
        Output.h)
