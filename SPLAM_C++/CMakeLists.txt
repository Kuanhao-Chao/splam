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
find_package(Python REQUIRED COMPONENTS Interpreter Development)

# find_package(Python COMPONENTS Interpreter Development)

# find_package (Python3 COMPONENTS Interpreter Development)

# find_package(PythonInterp 3.2 REQUIRED)
# find_package(PythonLibs 3.2 REQUIRED)

# set(PYTHON_INCLUDE_DIRS /Users/chaokuan-hao/miniconda3/include/python3.8/) # path to header
# set(PYTHON_LIBRARIES /Users/chaokuan-hao/miniconda3/lib/libpython3.8.dylib)

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
target_link_libraries(splam stdc++fs)
target_link_libraries(splam ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
# target_link_libraries(splam "${TORCH_LIBRARIES}")
# target_include_directories(splam PRIVATE ${PYTHON_INCLUDE_DIRS})

target_link_libraries(splam  Python::Python Python::Module)

# target_link_libraries(splam ${PYTHON_LIBRARIES})
set_property(TARGET splam PROPERTY CXX_STANDARD 17)

# The following code block is suggested to be used on Windows.
# According to https://github.com/pytorch/pytorch/issues/25457,
# the DLLs need to be copied to avoid memory errors.
# if (MSVC)
#   file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
#   add_custom_command(TARGET splam
#                      POST_BUILD
#                      COMMAND ${CMAKE_COMMAND} -E copy_if_different
#                      ${TORCH_DLLS}
#                      $<TARGET_FILE_DIR:splam>)
# endif (MSVC)

message("${CMAKE_CURRENT_SOURCE_DIR}/build/: ")
message("${CMAKE_CURRENT_SOURCE_DIR}/splam/: ")

message("Printing ${CMAKE_INSTALL_PREFIX}: ")

message("CMAKE_CURRENT_SOURCE_DIR:${CMAKE_CURRENT_SOURCE_DIR}")
message("PYTHON_FOUND:${PYTHON_FOUND}")
message("PYTHON_VERSION:${PYTHON_VERSION}")
# message("Python_Development_FOUND:${Python_Development_FOUND}")
message("PYTHON_LIBRARIES:${PYTHON_LIBRARIES}")
message("PYTHON_INCLUDE_DIRS:${PYTHON_INCLUDE_DIRS}")
message("CMAKE_C_FLAGS:${CMAKE_C_FLAGS}")

# install(TARGETS splam DESTINATION ${CMAKE_INSTALL_PREFIX})

# install(TARGETS splam DESTINATION /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/splam)


# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/"
#         DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/splam/")

# install(DIRECTORY "myDir/"
#         DESTINATION "myDestination"
#         FILES_MATCHING PATTERN "*.h" )

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/build/splam ${CMAKE_CURRENT_SOURCE_DIR}/splam/splam COPYONLY)



