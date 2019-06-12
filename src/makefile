#Compiler and library linking for local system
COMPILER= gfortran

FLAGS = -O3

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)


all: transform

transform: $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ 

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

constants.o: constants.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

statistics.o: statistics.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

special.o: special.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

av18.o: av18.f90 types.o constants.o
	$(COMPILER) $(FLAGS) -c $<

pion_exchange.o: pion_exchange.f90 types.o constants.o special.o
	$(COMPILER) $(FLAGS) -c $<

c2m_transform.o: c2m_transform.f90 types.o special.o constants.o av18.o
	$(COMPILER) $(FLAGS) -c $<

c2m_basis.o: c2m_basis.f90 types.o constants.o
	$(COMPILER) $(FLAGS) -c $<

c2m_montecarlo.o: c2m_montecarlo.f90 types.o c2m_basis.o c2m_transform.o statistics.o
	$(COMPILER) $(FLAGS) -c $<

c2m_io.o: c2m_io.f90 types.o c2m_transform.o c2m_montecarlo.o constants.o
	$(COMPILER) $(FLAGS) -c $<

c2m_main.o: c2m_main.f90 types.o c2m_basis.o c2m_io.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod 

mrproper: clean 
	rm -rf transform
