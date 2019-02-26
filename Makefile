### complier name
COMP:=gfortran

## Put .o and .mod files here
OBJDIR=build

# Path to search for files
 
VPATH = $(OBJDIR)

## OPTIONS
OPTS = -ffree-line-length-none -fcheck=bounds -Wall -Wextra -fbacktrace  -ffpe-trap=invalid,zero,overflow -fcheck=all -J$(OBJDIR) -I$(OBJDIR)

### All object files are here. If you want to add new .o files add to OBJ
OBJ = $(addprefix $(OBJDIR)/, Read_init.o Inputs.o Supermodel_main.o)

EXE = main

$(OBJDIR)/%.o: %.f90
	$(COMP) $(OPTS) -c $< -o $@

$(EXE): $(OBJ)
	$(COMP) $(OPTS) $(OBJ) -o $@


clean: 
	@cd $(OBJDIR) ; rm *
	rm -f *.o *.mod $(EXE)


