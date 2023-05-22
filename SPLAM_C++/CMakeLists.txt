cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(splam)
set(CMAKE_CXX_STANDARD 17)

# list(APPEND CMAKE_PREFIX_PATH "python -c 'import torch;print(torch.utils.cmake_prefix_path)'")
# find_package(Torch REQUIRED)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -fpermissive -DNOCURL=1")
# set(CMAKE_OSX_ARCHITECTURES "arm64")  # https://github.com/pytorch/pytorch/issues/63558

# set(Python_ADDITIONAL_VERSIONS 3.4)


# set(PYTHON_INCLUDE_DIR /Users/chaokuan-hao/miniconda3/include/python3.8/) # path to header
# set(PYTHON_LIBRARY /Users/chaokuan-hao/miniconda3/lib)

# find_package(PythonLibs REQUIRED)
# find_package(Python COMPONENTS Interpreter Development)
find_package(Python REQUIRED COMPONENTS Interpreter Development REQUIRED)

# find_package(Python COMPONENTS Interpreter Development)
# find_package (Python3 COMPONENTS Interpreter Development)
# find_package(PythonInterp 3.2 REQUIRED)
# find_package(PythonLibs 3.2 REQUIRED)

if(UNIX AND APPLE)
        set(PYTHON_INCLUDE_DIRS /Users/chaokuan-hao/miniconda3/include/python3.8/) # path to header
        set(PYTHON_LIBRARIES /Users/chaokuan-hao/miniconda3/lib/libpython3.8.dylib)
endif()

find_package(ZLIB)


# export LDFLAGS="-L/usr/local/opt/zlib/lib -L/usr/local/opt/bzip2/lib/"
# export CPPFLAGS="-I/usr/local/opt/zlib/include"
# export CPPFLAGS="-I/usr/local/opt/zlib/include"

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/zlib/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/local/opt/bzip2/include")




include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")
include_directories(${PYTHON_INCLUDE_DIRS})


include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoconf -i COMMAND ./configure  --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(splam
        src/splam.cpp
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

add_dependencies(splam htslib)
if(UNIX AND NOT APPLE)
        target_link_libraries(splam stdc++fs)
endif ()
target_link_libraries(splam ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# target_link_libraries(splam "${TORCH_LIBRARIES}")
# target_includlse_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})

if(UNIX AND NOT APPLE)
        target_link_libraries(splam  Python::Python Python::Module)
endif()

if(UNIX AND APPLE)
        target_link_libraries(splam ${PYTHON_LIBRARIES})
endif()

set_property(TARGET splam PROPERTY CXX_STANDARD 17)

message("${CMAKE_CURRENT_SOURCE_DIR}/build/: ")
message("${CMAKE_CURRENT_SOURCE_DIR}/splam/: ")

message("Printing ${CMAKE_INSTALL_PREFIX}: ")

message("CMAKE_CURRENT_SOURCE_DIR:${CMAKE_CURRENT_SOURCE_DIR}")
message("PYTHON_FOUND:${PYTHON_FOUND}")
message("Python_Development_FOUND:${Python_Development_FOUND}")
message("PYTHON_LIBRARIES:${PYTHON_LIBRARIES}")
message("PYTHON_INCLUDE_DIRS:${PYTHON_INCLUDE_DIRS}")
message("CMAKE_C_FLAGS:${CMAKE_C_FLAGS}")
