#Compiler and library linking for local system
COMPILER= gfortran

FLAGS = 

SRC = types.f90 constants.f90 special.f90 c2m_transform.f90 c2m_io.f90 c2m_basis.f90 c2m_main.f90

OBJ = $(SRC:.f90=.o)


all: transform

transform: $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ 

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

constants.o: constants.f90 types.mod
	$(COMPILER) $(FLAGS) -c $<

special.o: special.f90 types.mod
	$(COMPILER) $(FLAGS) -c $<

c2m_transform.o: c2m_transform.f90 types.mod special.mod
	$(COMPILER) $(FLAGS) -c $<

c2m_io.o: c2m_io.f90 types.mod c2m_transform.mod
	$(COMPILER) $(FLAGS) -c $<

c2m_basis.o: c2m_basis.f90 types.mod 
	$(COMPILER) $(FLAGS) -c $<

c2m_main.o: c2m_main.f90 types.mod c2m_io.mod c2m_basis.mod
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod 

mrproper: clean 
	rm -rf transform
