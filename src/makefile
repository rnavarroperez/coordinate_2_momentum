#Compiler and library linking for local system
COMPILER= gfortran

FLAGS = 

SRC = types.f90 c2m_parameters.f90 c2m_io.f90 c2m_main.f90

OBJ = $(SRC:.f90=.o)


all: transform

transform: $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ 

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

c2m_parameters.o: c2m_parameters.f90 types.mod
	$(COMPILER) $(FLAGS) -c $<

c2m_io.o: c2m_io.f90 types.mod c2m_parameters.mod
	$(COMPILER) $(FLAGS) -c $<

c2m_main.o: c2m_main.f90 types.mod c2m_parameters.mod c2m_io.mod
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod 

mrproper: clean 
	rm -rf transform
