FC = caf
FFLAGS = -Wall -O0 -g -fcheck=all
#FFLAGS = -Wall -O3 -g

progs = $(basename $(wildcard *.f90))

all: $(progs)

%: %.f90 Makefile
	$(FC) $(FFLAGS) -o $@ ~/lammps-velocityaverage/readers/reader.f90 $<

clean:
	rm -rf *.dat $(progs)
