cmake_minimum_required(VERSION 3.22)

project(
 DifSolver
 VERSION 0.01
 LANGUAGES C CXX)

add_compile_options(-O3)

add_subdirectory(sources)

set(include_dir ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(sources PRIVATE ${include_dir})

add_executable(main_scan main_scan.cpp)
add_executable(main_grad main_grad.cpp)

target_link_libraries(main_scan PRIVATE sources)
target_link_libraries(main_grad PRIVATE sources)