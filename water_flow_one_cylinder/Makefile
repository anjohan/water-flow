FC = gfortran
FFLAGS = -Wall -O3

all: vx_analysis

%: %.o
	$(FC) $(FFLAGS) -o $@ $^
%.o: %.f90
	$(FC) $(FFLAGS) -c $^
