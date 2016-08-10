FORT = gfortran
FLAGS = -march=native -ffast-math -funroll-loops -fstrict-aliasing -O3 -Wunused -cpp
SRCS = mGlobal.f90 mRandom.f90 mConfig.f90 EmDee.f90
OBJS = $(SRCS:%.f90=%.o)
LIBFILE = libmscommon.a
NSY=1

all: lib

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	ar -cr $(LIBFILE) $(OBJS)

EmDee.o: EmDee.f90
	$(FORT) $(FLAGS) -c -o $@ $<

mRandom.o: mRandom.f90
	$(FORT) $(FLAGS) -c -o $@ $<

mConfig.o: mConfig.f90 mGlobal.o
	$(FORT) $(FLAGS) -c -o $@ $<

mGlobal.o: mGlobal.f90
	$(FORT) $(FLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(LIBFILE)
	rm -rf html/

