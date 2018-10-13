cmake_minimum_required(VERSION 3.9)
project(virtualscannerexe VERSION 18.10.09 LANGUAGES CXX)

add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

add_subdirectory(../libvirtualscanner "${CMAKE_CURRENT_BINARY_DIR}/virtual_scanner")
add_executable(virtualscannerexe src/main.cpp)

target_compile_features(virtualscannerexe PRIVATE cxx_auto_type)

find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(OpenMP REQUIRED)

set_target_properties(virtualscannerexe PROPERTIES OUTPUT_NAME "virtualscanner")
target_compile_options(virtualscannerexe
    PRIVATE 
        $<$<CXX_COMPILER_ID:GNU>:-Wall -fPIC -Wno-sign-compare -Wno-uninitialized -frounding-math>
        $<$<CXX_COMPILER_ID:MSVC>: /Wall>)

target_link_libraries(virtualscannerexe
    PRIVATE
        OpenMP::OpenMP_CXX
        Boost::filesystem
        VirtualScanner::VirtualScanner)

