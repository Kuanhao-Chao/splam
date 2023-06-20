cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(splam)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(pybind11)


# pybind11_add_module(bind_test src/bind_test.cpp)
# target_compile_definitions(bind_test
#                            PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


find_package(ZLIB)

# export LDFLAGS="-L/usr/local/opt/zlib/lib -L/usr/local/opt/bzip2/lib/"
# export CPPFLAGS="-I/usr/local/opt/zlib/include"
# export CPPFLAGS="-I/usr/local/opt/zlib/include"

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/zlib/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/bzip2/include")




include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")
# include_directories(${PYTHON_INCLUDE_DIRS})


include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoconf -i COMMAND ./configure  --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )


############################
# pybind11 for building splam_extract 
############################
pybind11_add_module(splam_extract 
        src/splam_extract.cpp
        src/GSam.h
        src/extract.cpp
        src/junc_func.cpp
        src/util.cpp
        src/tmerge.cpp
        src/GSam.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)
target_link_libraries(splam_extract PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# set_property(TARGET splam_extract PROPERTY CXX_STANDARD 17)
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(splam_extract
           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


############################
# pybind11 for building splam-clean
############################
# Compiling splam-clean
pybind11_add_module(splam_clean
        src/splam_clean.cpp
        src/bam_sort.c
        src/bam.c
        src/sam_utils.c
        src/sam_opts.c
        src/bedidx.c
        src/junc_func.cpp
        src/extract.cpp
        src/predict.cpp
        src/clean.cpp
        src/update.cpp
        src/util.cpp
        src/tmerge.cpp
        src/GSam.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)
# target_link_libraries(splam Python::Python)
# add_dependencies(splam-clean htslib)
# if(UNIX AND NOT APPLE)
#         target_link_libraries(splam-clean stdc++fs)
# endif ()
target_link_libraries(splam_clean PRIVATE ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
target_compile_definitions(splam_clean
           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})






# ############################
# # Target for building splam-extract
# ############################
# # Compiling splam-extract
# add_executable(splam-extract
#         src/splam_extract.cpp
#         src/extract.cpp
#         src/junc_func.cpp
#         src/util.cpp
#         src/tmerge.cpp
#         src/GSam.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
#         ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

# # target_link_libraries(splam Python::Python)

# add_dependencies(splam-extract htslib)
# if(UNIX AND NOT APPLE)
#         target_link_libraries(splam-extract stdc++fs)
# endif ()
# target_link_libraries(splam-extract ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# # target_link_libraries(splam "${TORCH_LIBRARIES}")
# # target_includlse_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})

# set_property(TARGET splam-extract PROPERTY CXX_STANDARD 17)







############################
# Target for building splam-test
############################
# Compiling splam-test
add_executable(splam-test
        src/splam_test.cpp
        src/bam_sort.c
        src/bam.c
        src/sam_utils.c
        src/sam_opts.c
        src/bedidx.c)


# target_link_libraries(splam Python::Python)

add_dependencies(splam-test htslib)
if(UNIX AND NOT APPLE)
        target_link_libraries(splam-test stdc++fs)
endif ()
target_link_libraries(splam-test ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# target_link_libraries(splam "${TORCH_LIBRARIES}")
# target_includlse_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})


set_property(TARGET splam-test PROPERTY CXX_STANDARD 14)





# ############################
# # Target for building splam-merge
# ############################
# # Compiling splam-merge
# add_executable(splam-merge
#         src/bam_sort.c
#         src/bam.c
#         src/sam_utils.c
#         src/sam_opts.c
#         src/bedidx.c)

# # target_link_libraries(splam Python::Python)

# add_dependencies(splam-merge htslib)
# if(UNIX AND NOT APPLE)
#         target_link_libraries(splam-merge stdc++fs)
# endif ()
# target_link_libraries(splam-merge ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# # target_link_libraries(splam "${TORCH_LIBRARIES}")
# # target_includlse_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})

# if(UNIX AND NOT APPLE)
#         target_link_libraries(splam-merge  Python::Python Python::Module)
# endif()

# if(UNIX AND APPLE)
#         target_link_libraries(splam-merge ${PYTHON_LIBRARIES})
# endif()

# set_property(TARGET splam-merge PROPERTY CXX_STANDARD 17)













# ############################
# # pybind11 build splam-extract
# ############################
# # Compiling splam-extract
# add_executable(splam-extract
#         src/splam_extract.cpp
#         src/extract.cpp
#         src/junc_func.cpp
#         src/util.cpp
#         src/tmerge.cpp
#         src/GSam.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GHash.hh
#         ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
#         ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
#         ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

# # target_link_libraries(splam Python::Python)

# add_dependencies(splam-extract htslib)
# if(UNIX AND NOT APPLE)
#         target_link_libraries(splam-extract stdc++fs)
# endif ()
# target_link_libraries(splam-extract ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# # target_link_libraries(splam "${TORCH_LIBRARIES}")
# # target_includlse_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})

# set_property(TARGET splam-extract PROPERTY CXX_STANDARD 17)


# add_subdirectory(pybind11)
# pybind11_add_module(bind_test src/bind_test.cpp)

# target_compile_definitions(bind_test
#                            PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


message("${CMAKE_CURRENT_SOURCE_DIR}/build/: ")
message("${CMAKE_CURRENT_SOURCE_DIR}/splam/: ")
message("Printing ${CMAKE_INSTALL_PREFIX}: ")

message("CMAKE_CURRENT_SOURCE_DIR:${CMAKE_CURRENT_SOURCE_DIR}")
message("CMAKE_C_FLAGS:${CMAKE_C_FLAGS}")
