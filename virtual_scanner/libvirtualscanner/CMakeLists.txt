cmake_minimum_required(VERSION 3.9)
project(libvirtualscanner VERSION 18.10.09 LANGUAGES CXX)

if(NOT VCPKG_TOOLCHAIN)
    list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
endif()

##################################################
# Dependencies
find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
##################################################

add_library(virtualscanner
    src/virtual_scanner.cpp
    src/points.cpp)

add_library(VirtualScanner::VirtualScanner ALIAS virtualscanner)

target_include_directories(virtualscanner
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_features(virtualscanner PRIVATE cxx_nullptr)
target_compile_options(virtualscanner
    PRIVATE 
        $<$<CXX_COMPILER_ID:GNU>:-Wall -fPIC -Wno-sign-compare -Wno-uninitialized -frounding-math>
        $<$<CXX_COMPILER_ID:MSVC>: /Wall>)

target_link_libraries(virtualscanner
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        Boost::system
        CGAL::CGAL)

##################################################
# Installation

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/VirtualScanner)

install(TARGETS virtualscanner
    EXPORT virtualscanner-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties(virtualscanner PROPERTIES EXPORT_NAME VirtualScanner)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT virtualscanner-targets
    FILE
      VirtualScannerTargets.cmake
    NAMESPACE
      VirtualScanner::
    DESTINATION
      ${INSTALL_CONFIGDIR}
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/VirtualScannerConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/VirtualScannerConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/VirtualScannerConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
set(config_files)
list(APPEND config_files
    ${CMAKE_CURRENT_BINARY_DIR}/VirtualScannerConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/VirtualScannerConfigVersion.cmake
)
if(NOT VCPKG_TOOLCHAIN)
    list(APPEND config_files
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCGAL.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindEigen3.cmake
    )
endif()

install(FILES
    config_files
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree
if (NOT VCPKG_TOOLCHAIN)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindEigen3.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindEigen3.cmake
        COPYONLY)

    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCGAL.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindCGAL.cmake
        COPYONLY)
endif()

export(EXPORT virtualscanner-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/VirtualScannerTargets.cmake NAMESPACE VirtualScanner::)

#Register package in user's package registry
export(PACKAGE VirtualScanner)
