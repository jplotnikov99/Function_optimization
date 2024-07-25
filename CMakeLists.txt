cmake_minimum_required(VERSION 3.22)

project(
 DifSolver
 VERSION 0.01
 LANGUAGES C CXX)

add_compile_options(-O3)

add_subdirectory(sources)
add_subdirectory(functions_sources)

set(include_dir1 ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(include_dir2 ${CMAKE_CURRENT_SOURCE_DIR}/functions)

target_include_directories(sources PRIVATE ${include_dir1} ${include_dir2})
target_include_directories(functions_sources PRIVATE ${include_dir1} ${include_dir2})

add_executable(main_scan main_scan.cpp)
add_executable(main_grad main_grad.cpp)

target_link_libraries(main_scan PRIVATE sources functions_sources)
target_link_libraries(main_grad PRIVATE sources functions_sources)