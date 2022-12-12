cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(intron_matcher)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DNOCURL=1")

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

add_executable(intron_matcher src/intron_matcher.cpp 
        src/GSam.cpp
        src/tmerge.cpp
        src/commons.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.h

        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.h

        ${CMAKE_SOURCE_DIR}/include/gclib/codons.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GList.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GHashMap.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GBitVec.h
        
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

# add_dependencies(intron_matcher htslib)
target_link_libraries(intron_matcher ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
        
set_property(TARGET intron_matcher PROPERTY CXX_STANDARD 17)


install(TARGETS intron_matcher DESTINATION ${CMAKE_INSTALL_PREFIX})