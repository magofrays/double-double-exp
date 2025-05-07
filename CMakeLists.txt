cmake_minimum_required(VERSION 3.25)
project(DoubleDouble LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  

add_executable(
    DoubleDouble
    src/headers/doubledouble.cpp
    src/main.cpp   
)

include_directories(DoubleDouble PUBLIC src/headers)