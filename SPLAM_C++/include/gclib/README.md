## GCLib - Genomic C++ Library
This is an eclectic collection of basic C++ code (functions, classes, templates) which is shared between a few of my bioinformatics projects. The main idea was to provide a core collection of data structures, trying to avoid unnecessary code dependencies of other heavy libraries, while minimizing build time. 

I had started gathering this code even before the C++ STL had been fully adopted as a cross-platform "standard". Even STL itself seems a bit on the heavy side (and keeps growing) compared to what I need in practice for many of my C++ projects, so often times I prefer to just use these simpler and leaner C++ classes and templates to provide most common data structures needed for my projects.

## Build/Install
Do not build. Do not install. This is not meant to be built into an object library, it's a simple _source code library_ for other projects to include and link statically into the final executable(s). The makefile included here is just for simple, extemporaneous tests I occasionally perform as new functionality is added to this code collection.

