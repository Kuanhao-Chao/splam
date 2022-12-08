cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SpliceNN)
set(CMAKE_CXX_STANDARD 17)

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -fpermissive -DNOCURL=1")

include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ./configure --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(SpliceNN src/SpliceNN.cpp 
        src/junc.cpp 
        src/GSam.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

add_dependencies(SpliceNN htslib)
target_link_libraries(SpliceNN ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
        
target_link_libraries(SpliceNN "${TORCH_LIBRARIES}")
set_property(TARGET SpliceNN PROPERTY CXX_STANDARD 17)


# The following code block is suggested to be used on Windows.
# According to https://github.com/pytorch/pytorch/issues/25457,
# the DLLs need to be copied to avoid memory errors.
if (MSVC)
  file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
  add_custom_command(TARGET SpliceNN
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${TORCH_DLLS}
                     $<TARGET_FILE_DIR:SpliceNN>)
endif (MSVC)

install(TARGETS SpliceNN DESTINATION ${CMAKE_INSTALL_PREFIX})