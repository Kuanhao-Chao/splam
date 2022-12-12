# the compiler: gcc for C program, define as g++ for C++
CXX = g++
# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS = -std=c++17 -Wall -O3 -pthread
GDIR := ./gclib
INCDIRS := -I. -I${GDIR}

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GSam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

# The build target 
EXE = intron_matcher
 
$(EXE): intron_matcher.o $(OBJS)
	@echo "> Compiling $(EXE)"
	@$(CXX) $^ -o $(EXE) $(CFLAGS) $(INCLUDE) -w
	@cp $(EXE) ../../../bin
	@echo $(EXE) copied to ../../../bin

%.o: %.cpp
	@echo "> Compiling $^"
	@$(CXX) $(CFLAGS) -c $(INCLUDE)  $<

remove:
	@rm -R -- */

clean:
	@rm -R -- */ $(TARGET) *.o