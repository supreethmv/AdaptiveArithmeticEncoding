GPP=g++
CCFLAGS=-msse -Wall
CCFLAGS+=-g
#CCFLAGS+=-O3
LDFLAGS=-lm 

OBJECTS=src/bfio.o \
	src/image_io.o \
	src/alloc.o 

all: compress

compress: $(OBJECTS) src/ic18_adaptive_arithmetic.c Makefile
	$(GPP) $(CCFLAGS) $(OBJECTS) src/ic18_adaptive_arithmetic.c -o ic18_adaptive_arithmetic $(LDFLAGS) 

%.o : %.c
	$(GPP) $(CCFLAGS) -o $@ -c $<

clean:
	rm -f ic18_adaptive_arithmetic
	rm -f src/*.o src/*~
	rm -f *.o *~
