cmake_minimum_required(VERSION 3.15...3.26)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(nautilus LANGUAGES CXX)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++")

if (NOT SKBUILD)
    message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve -Ceditable.rebuild=true .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

add_subdirectory(fftw)
add_subdirectory(rfftw)
add_subdirectory(ccp4)
add_subdirectory(mmdb2)
add_subdirectory(gemmi)
add_subdirectory(clipper/gemmi)
add_subdirectory(clipper/minimol)
add_subdirectory(clipper/core)
add_subdirectory(clipper/contrib)
add_subdirectory(clipper/cns)
add_subdirectory(clipper/mmdb)
add_subdirectory(clipper/cif)
add_subdirectory(clipper/phs)
add_subdirectory(clipper/ccp4)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        checkout/mmdb2
        checkout/clipper
        checkout/libccp4/ccp4
        checkout/gemmi/include
)

link_directories(
        mmdb2
        clipper
        ccp4
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)


add_library(
        nautilus_lib
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-join.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-mlfind.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-rebuild-bases.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-refine.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-sequence.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-ss-find.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-target.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-tidy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-tools.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-util.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nucleicacid_db.cpp
)


add_executable(cnautilus ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/cnautilus.cpp)
set_target_properties(cnautilus PROPERTIES OUTPUT_NAME nautilus)
target_compile_features(cnautilus INTERFACE cxx_std_17)

target_include_directories(cnautilus PUBLIC checkout/gemmi)
target_link_directories(cnautilus
        PUBLIC
        mmdb2
        clipper/mmdb
        clipper/cif
        clipper/core
        clipper/contrib
        clipper/ccp4
        clipper/gemmi
        clipper/minimol
        fftw
        rfftw
        ccp4
)


target_link_libraries(
        cnautilus
        PRIVATE
        gemmi_cpp
        nautilus_lib
        ccp4
        fftw
        rfftw
        mmdb2
        clipper-mmdb
        clipper-cif
        clipper-core
        clipper-contrib
        clipper-minimol
        clipper-ccp4
        clipper-gemmi
)
set(nautilus_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-join.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-mlfind.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-rebuild-bases.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-refine.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-sequence.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-ss-find.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-target.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-tidy.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-tools.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nautilus-util.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/nucleicacid_db.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/python/nautilus-include.h

)
set_target_properties(cnautilus PROPERTIES PUBLIC_HEADER "${nautilus_HEADERS}")
install(TARGETS cnautilus
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/nautilus)


# Try to import all Python components potentially needed by nanobind
find_package(Python 3.8
        REQUIRED COMPONENTS Interpreter Development.Module
        OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# We are now ready to compile the actual extension module
nanobind_add_module(
        nautilus_module

        STABLE_ABI
        NB_STATIC

        # Python bound functions here
        src/python/pynautilus.cpp
)

target_include_directories(nautilus_module PUBLIC
        checkout/gemmi
        mmdb2
        clipper/mmdb
        clipper/cif
        clipper/core
        clipper/contrib
        clipper/ccp4
        clipper/gemmi
        clipper/minimol
        fftw
        rfftw
        ccp4
)


target_link_directories(nautilus_module
        PUBLIC
        mmdb2
        clipper/mmdb
        clipper/cif
        clipper/core
        clipper/contrib
        clipper/ccp4
        clipper/gemmi
        clipper/minimol
        fftw
        rfftw
        ccp4
)


target_link_libraries(
        nautilus_module
        PRIVATE
        gemmi_cpp
        nautilus_lib
        ccp4
        fftw
        rfftw
        mmdb2
        clipper-mmdb
        clipper-cif
        clipper-core
        clipper-contrib
        clipper-minimol
        clipper-ccp4
        clipper-gemmi
)
# Install directive for scikit-build-core
install(TARGETS nautilus_module LIBRARY DESTINATION pynautilus)
