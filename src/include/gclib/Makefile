INCDIRS := 
#-I${BAM}

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)
LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)
LIBS    := -lz

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (, $(findstring mingw32, $(DMACH)))
WINDOWS=1
endif

ifneq (, $(findstring linux, $(DMACH)))
 # -lrt only needed for Linux systems
 LIBS+= -lrt
endif


# MinGW32 GCC 4.5 link problem fix
ifdef WINDOWS
 DMACH := windows
 ifeq ($(findstring 4.5.,$(shell ${CXX} -dumpversion)), 4.5.)
  LDFLAGS += -static-libstdc++ -static-libgcc
 endif
  LIBS += -lregex -lws2_32
endif

# Misc. system commands
# assuming on Windows this is run under MSYS
RM = rm -f

# File endings
ifdef WINDOWS
EXE = .exe
else
EXE =
endif

BASEFLAGS  := -std=c++11 -Wall -Wextra ${INCDIRS} $(MARCH) \
 -D_REENTRANT -fno-exceptions -fno-rtti

GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

#add the link-time optimization flag if gcc version > 4.5
GCC_VERSION:=$(subst ., ,$(shell gcc -dumpversion))
GCC_MAJOR:=$(word 1,$(GCC_VERSION))
GCC_MINOR:=$(word 2,$(GCC_VERSION))
#GCC_SUB:=$(word 3,$(GCC_VERSION))
GCC_SUB:=x

ifneq (,$(filter %release %nodebug, $(MAKECMDGOALS)))
  # -- release build
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %tsan, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
   else
   # -- plain debug build
    CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS), -O0 -g)
    ifneq (, $(findstring darwin, $(DMACH)))
        CXXFLAGS += -gdwarf-3
    endif
    CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
   endif
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

.PHONY : all
all:    htest
memcheck tsan: all
#mdtest
nodebug: all
release: all
debug: all

OBJS := GBase.o GStr.o GArgs.o GResUsage.o

version: ; @echo "GCC Version is: "$(GCC_MAJOR)":"$(GCC_MINOR)":"$(GCC_SUB)
htest.o: htest.cpp GHashMap.hh
htest:  $(OBJS) htest.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
mdtest: $(OBJS) mdtest.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
# target for removing all object files

.PHONY : clean
clean:: 
	${RM} $(OBJS) *.o mdtest$(EXE) htest$(EXE)
	${RM} core.*
