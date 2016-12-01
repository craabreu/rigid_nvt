FORT = gfortran
FLAGS = -fno-range-check -march=native -ffast-math -funroll-loops -fstrict-aliasing -O3 -Wunused -cpp -fopenmp
COMMON = common
SRCS = $(wildcard $(COMMON)/*.f90)
EXECUTABLE ?= expens
LIB = mscommon
LIBFILE = $(COMMON)/lib$(LIB).a

.PHONY: all lj lj_coul

all: lj_coul

lj: $(EXECUTABLE)_lj

lj_coul: $(EXECUTABLE)_lj_coul

$(EXECUTABLE)_lj_coul: $(EXECUTABLE).f90 $(LIBFILE) $(SRCS)
	$(FORT) $(FLAGS) -Dcoul -J$(COMMON) -L$(COMMON) -o $@ $< -l$(LIB) -lemdee

$(EXECUTABLE)_lj: $(EXECUTABLE).f90 $(LIBFILE) $(SRCS)
	$(FORT) $(FLAGS) -J$(COMMON) -L$(COMMON) -o $@ $< -l$(LIB) -lemdee

$(LIBFILE): $(SRCS)
	make -C $(COMMON) lib

clean:
	rm -f $(EXECUTABLE)_lj
	rm -f $(EXECUTABLE)_lj_coul

clean-all:
	rm -f $(EXECUTABLE)_lj
	rm -f $(EXECUTABLE)_lj_coul
	make -C $(COMMON) clean

