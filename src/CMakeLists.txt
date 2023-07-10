cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(splam)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(pybind11)

find_package(ZLIB)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/zlib/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/bzip2/include")

include_directories("${CMAKE_SOURCE_DIR}/include/")
# link_directories("${CMAKE_SOURCE_DIR}/include/htslib")
message("${CMAKE_SOURCE_DIR}/include/htslib")

message("LIBRARY_OUTPUT_DIRECTORY: ${LIBRARY_OUTPUT_DIRECTORY}")


include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        # CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
        CONFIGURE_COMMAND autoreconf -i
        COMMAND ./configure --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND make CFLAGS=-fPIC
        INSTALL_COMMAND ""
        )

############################
# pybind11 for building splam_extract
############################
pybind11_add_module(splam_extract
        splam/module/splam_extract.cpp
        splam/module/GSam.h
        splam/module/extract.cpp
        splam/module/junc_func.cpp
        splam/module/util.cpp
        splam/module/tmerge.cpp
        splam/module/GSam.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)



# Add include directories
target_include_directories(splam_extract PRIVATE ${CMAKE_SOURCE_DIR}/include/)

# Link the library to your program
add_dependencies(splam_extract htslib)
if(UNIX AND NOT APPLE)
        target_link_libraries(splam_extract PRIVATE stdc++fs)
#         target_link_libraries(splam_extract PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# elseif (APPLE)
#         target_link_libraries(splam_extract PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
endif ()
target_link_libraries(splam_extract PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)


# define (VERSION_INFO) here.
target_compile_definitions(splam_extract
           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


############################
# pybind11 for building splam-clean
############################
pybind11_add_module(splam_clean
        splam/module/splam_clean.cpp
        splam/module/bam_sort.c
        splam/module/bam.c
        splam/module/sam_utils.c
        splam/module/sam_opts.c
        splam/module/bedidx.c
        splam/module/junc_func.cpp
        splam/module/extract.cpp
        splam/module/clean.cpp
        splam/module/update.cpp
        splam/module/util.cpp
        splam/module/tmerge.cpp
        splam/module/GSam.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

# Add include directories
target_include_directories(splam_clean PRIVATE ${CMAKE_SOURCE_DIR}/include/)

# Link the library to your program
add_dependencies(splam_clean htslib)
if(UNIX AND NOT APPLE)
        target_link_libraries(splam_clean PRIVATE stdc++fs)
# elseif (APPLE)
endif ()
target_link_libraries(splam_clean PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)


target_compile_definitions(splam_clean
           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

message("EXAMPLE_VERSION_INFO: ${EXAMPLE_VERSION_INFO}")


# # Move the build ${CMAKE_SOURCE_DIR}/include/htslib/libhts.so / ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a
# if(UNIX AND NOT APPLE)
#         configure_file(${CMAKE_SOURCE_DIR}/include/htslib/libhts.so ${CMAKE_SOURCE_DIR}/include/htslib/libhts.so COPYONLY)
# elseif (APPLE)
#         configure_file(${CMAKE_SOURCE_DIR}/include/htslib/libhts.a ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a COPYONLY)
# endif ()
