project( Code_ )

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)

FIND_PACKAGE(OpenMP)

include_directories(/home/mkerber/phat_swap_reduction/include)

find_package(Boost REQUIRED timer)

set (CMAKE_CXX_FLAGS "-fopenmp")

add_executable(main "main.cpp")
target_link_libraries(main "-lboost_timer")
