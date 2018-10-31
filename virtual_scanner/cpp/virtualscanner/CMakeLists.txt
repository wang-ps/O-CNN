cmake_minimum_required(VERSION 3.9)

add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

add_executable(virtualscannerexe src/main.cpp)

target_compile_features(virtualscannerexe PRIVATE cxx_auto_type)

find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(OpenMP REQUIRED)

set_target_properties(virtualscannerexe PROPERTIES
      OUTPUT_NAME "virtualscanner"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

target_compile_options(virtualscannerexe
    PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -fPIC -Wno-sign-compare -Wno-uninitialized -frounding-math>
        $<$<CXX_COMPILER_ID:MSVC>: /Wall>)

target_link_libraries(virtualscannerexe
    PRIVATE
        OpenMP::OpenMP_CXX
        Boost::filesystem
        OCNN::VirtualScanner)

