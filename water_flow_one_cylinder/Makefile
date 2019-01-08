FC = caf
#FFLAGS = -Wall -O0 -g -fcheck=all
FFLAGS = -Wall -O3 -g

progs = $(basename $(wildcard *.f90))

all: $(progs)

%: %.f90 Makefile
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -rf *.dat $(progs)
