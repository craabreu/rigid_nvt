FORT = gfortran
FLAGS = -march=native -ffast-math -funroll-loops -fstrict-aliasing -O3 -Wunused -cpp
SRCS = $(wildcard *.f90)
OBJS = $(SRCS:%.f90=%.o)
LIBFILE = libmscommon.a

all: lib

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	ar -cr $(LIBFILE) $(OBJS)

%.o: %.f90
	$(FORT) $(FLAGS) -c -o $@ $<

mThermostat.o: mThermostat.f90 mGlobal.o
	$(FORT) $(FLAGS) -c -o $@ $<

mGlobal.o: mGlobal.f90
	$(FORT) $(FLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(LIBFILE)
	rm -rf html/

