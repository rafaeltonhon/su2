FC=gfortran
FFLAGS=-O3 -Wall -fcheck=bounds

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

OBJECTS=\
   su2main.o\
   su2lattice.o\
   su2thermalize.o\
   su2measures.o\
   su2center.o\

SRC=\
   su2main.f90\
   su2lattice.f90\
   su2thermalize.f90\
   su2measures.f90\
   su2center.f90\

su2main: $(OBJECTS)
	$(FC) -o $@ $^

clean:
	rm *\.o su2main

ctags: $(SRC)
	ctags $(SRC)

su2main.o: su2lattice.o su2thermalize.o su2measures.o su2center.o su2main.f90
	$(FC) -Wall -c su2main.f90 -o su2main.o
	
su2center.o: su2lattice.o su2measures.o su2center.f90
	$(FC) -Wall -c su2center.f90 -o su2center.o
	
su2thermalize.o: su2lattice.o su2thermalize.f90
	$(FC) -Wall -c su2thermalize.f90 -o su2thermalize.o

measures.o: su2lattice.o su2measures.f90
	$(FC) -Wall -c su2measures.f90 -o su2measures.o
	
lattice.o: su2lattice.f90
	$(FC) -Wall -c su2lattice.f90 -o su2lattice.o
