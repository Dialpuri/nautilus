cmake_minimum_required(VERSION 3.12)
project(NAUTILUS LANGUAGES C CXX)

# set(CMAKE_VERBOSE_MAKEFILE TRUE)

SET(NAUTILUS_CC_FLAGS "-g -O0 -w")
SET(NAUTILUS_CXX_FLAGS "-g -O0 -w -fwasm-exceptions -I${CMAKE_INSTALL_PREFIX}/include --bind -lnodefs.js --profiling")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NAUTILUS_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NAUTILUS_CXX_FLAGS} ")
# add_compile_options(-std=c++11)

# option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "/Users/dialpuri/Development/privateer/emsdk/upstream/emscripten/em++")
add_link_options(-sDEFAULT_PTHREAD_STACK_SIZE=2MB -sSTACK_SIZE=5MB)
add_compile_options(--preload-file data@/)

add_subdirectory(fftw)
add_subdirectory(rfftw)
add_subdirectory(ccp4)
add_subdirectory(mmdb2)
add_subdirectory(clipper/minimol)
add_subdirectory(clipper/core)
add_subdirectory(clipper/ccp4)
add_subdirectory(clipper/contrib)
add_subdirectory(clipper/cns)
add_subdirectory(clipper/mmdb)
add_subdirectory(clipper/cif)
add_subdirectory(clipper/phs)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    checkout/mmdb2 
    checkout/clipper 
)

link_directories(
    checkout/mmdb2
    checkout/clipper
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_library(
    nautilus_lib
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-join.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-rebuild-bases.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-sequence.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-ss-find.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-target.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tidy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tools.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-util.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nucleicacid_db.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-mlfind.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-refine.cpp
)


# add_executable(cnautilus
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/cnautilus.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-join.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-rebuild-bases.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-sequence.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-ss-find.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-target.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tidy.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tools.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-util.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nucleicacid_db.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-mlfind.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-refine.cpp
#     )

add_executable(cnautilus_exec ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-bind.cpp)
set_target_properties (cnautilus_exec PROPERTIES OUTPUT_NAME pnautilus )
target_link_options(cnautilus_exec PRIVATE --preload-file data@/ -sDEMANGLE_SUPPORT=1 -sASSERTIONS=2 -fexceptions -sALLOW_MEMORY_GROWTH=1 --bind -sFORCE_FILESYSTEM=1 -sMODULARIZE=1 -sEXPORT_NAME=nautilus_module -sEXPORTED_RUNTIME_METHODS=['FS'] -sMODULARIZE=1 -sEXPORT_ES6=1 --bind -lnodefs.js --profiling)

target_link_directories(cnautilus_exec
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
 )

target_link_libraries(
    cnautilus_exec
    nautilus_lib
    ccp4
    mmdb2
    clipper-mmdb
    clipper-cif
    clipper-core 
    clipper-contrib 
    clipper-ccp4 
    clipper-minimol
    fftw 
    rfftw 
    )

set(nautilus_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-join.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-rebuild-bases.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-sequence.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-ss-find.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-target.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tidy.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-tools.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-util.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nucleicacid_db.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-mlfind.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nautilus-refine.h
)

set_target_properties(cnautilus_exec PROPERTIES PUBLIC_HEADER "${nautilus_HEADERS}")
# set_target_properties(cnautilus_exec PROPERTIES COMPILE_FLAGS "-s ENVIRONMENT='web'")

install(TARGETS nautilus_lib
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/nautilus
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cnautilus.wasm
${CMAKE_CURRENT_BINARY_DIR}/cnautilus.js
DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/webserver/wasm)
