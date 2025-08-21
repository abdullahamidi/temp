cmake_minimum_required(VERSION 3.12)
project(cgns2foam VERSION 1.0.0 LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/.deps/.install")

# Add compiler flags for better debugging and warnings
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
add_compile_options(-march=native)

set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 

# Find required packages

find_package(Boost REQUIRED COMPONENTS log log_setup filesystem thread system program_options json)
find_package(TBB REQUIRED)

link_libraries(TBB::tbb)

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)

# Find CGNS library with flexible path support
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(CGNS QUIET cgns)
endif()

# If pkg-config didn't find CGNS, try manual search
if(NOT CGNS_FOUND)
    # Check for CGNS_ROOT environment variable or CMake variable
    if(DEFINED ENV{CGNS_ROOT})
        set(CGNS_ROOT $ENV{CGNS_ROOT})
    endif()
    
    if(CGNS_ROOT)
        set(CGNS_INCLUDE_HINTS ${CGNS_ROOT}/include)
        set(CGNS_LIBRARY_HINTS ${CGNS_ROOT}/lib ${CGNS_ROOT}/lib64)
    endif()
    
    # Find CGNS header
    find_path(CGNS_INCLUDE_DIR
        NAMES cgnslib.h
        HINTS ${CGNS_INCLUDE_HINTS}
        PATHS
            /usr/include
            /usr/local/include
            /opt/cgns/include
            /opt/local/include
    )
    
    # Find CGNS library
    find_library(CGNS_LIBRARY
        NAMES cgns CGNS
        HINTS ${CGNS_LIBRARY_HINTS}
        PATHS
            /usr/lib
            /usr/lib64
            /usr/local/lib
            /usr/local/lib64
            /opt/cgns/lib
            /opt/cgns/lib64
            /opt/local/lib
    )
    
    if(CGNS_INCLUDE_DIR AND CGNS_LIBRARY)
        set(CGNS_FOUND TRUE)
        set(CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIR})
        set(CGNS_LIBRARIES ${CGNS_LIBRARY})
    endif()
endif()

# Check if CGNS was found
if(NOT CGNS_FOUND)
    message(FATAL_ERROR "CGNS library not found. Please install CGNS or set CGNS_ROOT to the installation directory.")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CGNS_INCLUDE_DIRS})

# Create executable
file(GLOB CGNS2FOAM_SRC "src/*.cpp")

add_executable(cgns2foam ${CGNS2FOAM_SRC})



# Link libraries
target_link_libraries(cgns2foam 
    ${Boost_LIBRARIES}
    ${CGNS_LIBRARIES}
    ${HDF5_LIBRARIES}
)

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  Project: ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Boost Version: ${Boost_VERSION}")
message(STATUS "  Boost Libraries: ${Boost_LIBRARIES}")
message(STATUS "  CGNS Include: ${CGNS_INCLUDE_DIR}")
message(STATUS "  CGNS Library: ${CGNS_LIBRARIES}")
message(STATUS "  HDF5 Library: ${HDF5_LIBRARIES}")
message(STATUS "")


install(TARGETS cgns2foam DESTINATION bin)