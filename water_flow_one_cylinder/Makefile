FC = gfortran
FFLAGS = -Wall -O3 -fopenmp

all: vx_analysis vx_profile_analysis

%: %.o Makefile
	$(FC) $(FFLAGS) -o $@ $<
%.o: %.f90 Makefile
	$(FC) $(FFLAGS) -c $<
