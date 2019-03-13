# ------------------Macro-Defs---------------------

# compiler
F90 = gfortran

# Put .o and .mod files here:
OBJDIR = src/build

# When compiling, search for files in these directories:
VPATH = $(OBJDIR):src:src/ACDC/ACDC_module_2016_09_23:src/ACDC/ACDC_module_ions_2018_08_31

# Options reminders:
# -w suppresses warning messages

# For programming (faster compiling):
BOX_OPTS = -g -ffree-line-length-none -cpp -DLINUX -DNETOUT -J$(OBJDIR) -I$(OBJDIR) -fcheck=bounds,do -Wall -Wextra -Wsurprising \
-Warray-temporaries -Wtabs -Wno-tabs -Wno-character-truncation -fbacktrace -ffpe-trap=invalid,zero,overflow -pg -g -fcheck=all
#BOX_OPTS = -ffree-line-length-none -cpp -DLINUX -J$(OBJDIR) -I$(OBJDIR) -fcheck=bounds,do -Wall -Wextra -Wsurprising \
#-Warray-temporaries -Wtabs -Wno-character-truncation -fbacktrace -ffpe-trap=invalid,zero,overflow -pg -g -fcheck=all

ACDC_OPTS = -ffree-line-length-none -cpp -J$(OBJDIR) -I$(OBJDIR) -fcheck=all -ffpe-trap=invalid,zero,overflow -O3


BOX_OBJECTS = second_precision.o constants.o auxillaries.o 


ACDC_OBJECTS = get_acdc_J.o driver.o acdc_system.o acdc_equations.o monomer_settings.o solution_settings.o \
               vode.o vodea.o

ACDC_D_OBJECTS = get_acdc_D.o driverd.o acdc_systemd.o acdc_equationsD.o monomer_settingsd.o solution_settings.o \
                 vode.o vodea.o

NETLIBS =  -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdf  -lnetcdff -lcurl

all: superbox.exe

# Here is the link step:
superbox.exe: Superbox.o $(BOX_OBJECTS) $(ACDC_OBJECTS) $(ACDC_D_OBJECTS)
	$(F90) $(BOX_OPTS) $^ -o $@ $(NETLIBS)


# Here are the compile steps
# Main program

#$(OBJDIR)/Superbox.o: Superbox.f90 $(CHEM_OBJECTS) $(BOX_OBJECTS) $(UHMA_OBJECTS) $(MEGAN_OBJECTS)
$(OBJDIR)/Superbox.o: src/Supermodel_main.f90 $(BOX_OBJECTS) $(ACDC_OBJECTS) $(ACDC_D_OBJECTS)
	 $(F90) $(BOX_OPTS) -c $< -o $@



# # Chemistry
# $(OBJDIR)/second_Precision.o: chemistry/second_Precision.f90
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Parameters.o: chemistry/second_Parameters.f90 second_Precision.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Global.o: chemistry/second_Global.f90 second_Parameters.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Initialize.o: chemistry/second_Initialize.f90 second_Parameters.o second_Global.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Util.o: chemistry/second_Util.f90 second_Parameters.o second_Monitor.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Monitor.o: chemistry/second_Monitor.f90
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_JacobianSP.o: chemistry/second_JacobianSP.f90
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_LinearAlgebra.o: chemistry/second_LinearAlgebra.f90 second_Parameters.o second_JacobianSP.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Jacobian.o: chemistry/second_Jacobian.f90 second_Parameters.o second_JacobianSP.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Rates.o: chemistry/second_Rates.f90 second_Parameters.o second_Global.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Integrator.o: chemistry/second_Integrator.f90 second_Precision.o second_Global.o second_Parameters.o second_JacobianSP.o \
# second_LinearAlgebra.o second_Function.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Function.o: chemistry/second_Function.f90 second_Parameters.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Model.o: chemistry/second_Model.f90 second_Precision.o second_Parameters.o second_Global.o second_Function.o \
# second_Integrator.o second_Rates.o second_Jacobian.o second_LinearAlgebra.o second_Monitor.o second_Util.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@
# $(OBJDIR)/second_Main.o: chemistry/second_Main.f90 second_Model.o second_Initialize.o
# 	 $(F90) $(CHEM_OPTS) -c $< -o $@



# # Uhma (aerosols)
# $(OBJDIR)/uhma_datatypes.o: uhma/uhma_datatypes.f90
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/uhma_io.o: uhma/uhma_io.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/distribute.o: uhma/distribute.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/nucleate.o: uhma/nucleate.f90 uhma_datatypes.o $(ACDC_OBJECTS) $(ACDC_D_OBJECTS)
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/condense.o: uhma/condense.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/coagulate.o: uhma/coagulate.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/deposit.o: uhma/deposit.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/uhma_misc.o: uhma/uhma_misc.f90 uhma_datatypes.o arbitrary_input.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/gde_solver.o: uhma/gde_solver.f90 uhma_datatypes.o nucleate.o condense.o coagulate.o deposit.o distribute.o \
# uhma_misc.o sorting.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/output_netcdf.o: uhma/output_netcdf.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@ $(NETLIBS)
# $(OBJDIR)/sorting.o: uhma/sorting.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@
# $(OBJDIR)/arbitrary_input.o: uhma/arbitrary_input.f90 uhma_datatypes.o
# 	$(F90) $(UHMA_OPTS) -c $< -o $@

# ACDC
$(OBJDIR)/get_acdc_J.o: ACDC/ACDC_module_ions_2018_08_31/get_acdc_J_ions.f90 acdc_system.o driver.o vode.o vodea.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/driver.o: ACDC/ACDC_module_ions_2018_08_31/driver_acdc_J_ions.f90 acdc_system.o solution_settings.o \
                   monomer_settings.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/acdc_system.o: ACDC/ACDC_module_ions_2018_08_31/acdc_system_AN_ions.f90
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/acdc_equations.o: ACDC/ACDC_module_ions_2018_08_31/acdc_equations_AN_ions.f90 monomer_settings.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/monomer_settings.o: ACDC/ACDC_module_ions_2018_08_31/monomer_settings_acdc_NH3_ions.f90 acdc_system.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

# ACDC for DMA
$(OBJDIR)/get_acdc_D.o: ACDC/ACDC_module_2016_09_23/get_acdc_D.f90 acdc_systemd.o driverd.o vode.o vodea.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/driverd.o: ACDC/ACDC_module_2016_09_23/driver_acdc_D.f90 acdc_systemd.o solution_settings.o \
                   monomer_settingsd.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/acdc_systemd.o: ACDC/ACDC_module_2016_09_23/acdc_system_AD_new.f90
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/acdc_equationsD.o: ACDC/ACDC_module_2016_09_23/acdc_equations_AD_new.f90 monomer_settingsd.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/monomer_settingsd.o: ACDC/ACDC_module_2016_09_23/monomer_settings_acdc_DMA.f90 acdc_systemd.o
	 $(F90) $(ACDC_OPTS) -c $< -o $@

# Common solvers for ACDC
$(OBJDIR)/solution_settings.o: ACDC/ACDC_module_ions_2018_08_31/solvers/solution_settings.f90
	 $(F90) $(ACDC_OPTS) -c $< -o $@

$(OBJDIR)/vode.o: ACDC/ACDC_module_ions_2018_08_31/solvers/vode.f
	gfortran -std=legacy -O3 -c $< -o $@

$(OBJDIR)/vodea.o: ACDC/ACDC_module_ions_2018_08_31/solvers/vodea.f
	gfortran -std=legacy -O3 -c $< -o $@

# Actual model files
$(OBJDIR)/auxillaries.o: src/auxillaries.f90
	$(F90) $(BOX_OPTS) -c $< -o $@

$(OBJDIR)/constants.o: src/constants.f90
	$(F90) $(BOX_OPTS) -c $< -o $@

$(OBJDIR)/second_precision.o: src/second_precision.f90
	$(F90) $(BOX_OPTS) -c $< -o $@

BOX_MODS = $(BOX_OBJECTS:.o=.mod)

# UHMA_MODS = $(UHMA_OBJECTS:.o=.mod)
# #CHEM_MODS = $(CHEM_OBJECTS:.o=.mod)
# CHEM_MODS = second_precision.mod second_parameters.mod second_initialize.mod second_util.mod second_monitor.mod second_jacobiansp.mod \
#                second_linearalgebra.mod second_jacobian.mod second_global.mod second_rates.mod second_integrator.mod second_function.mod \
#                second_model.mod second_main.mod
ACDC_MODS = $(ACDC_OBJECTS:.o=.mod)
ACDC_D_MODS = $(ACDC_D_OBJECTS:.o=.mod)

# This entry allows you to type 'make clean' to get rid of all object and module files

# With 'clean', don't remove chemistry object files, since it takes very long (30 min.) to compile them,
# and there usually is no need to recompile them
clean:
	-@cd $(OBJDIR) ; rm $(BOX_OBJECTS) $(BOX_MODS)       2>/dev/null || true
	-@cd $(OBJDIR) ; rm $(ACDC_OBJECTS) $(ACDC_MODS)     2>/dev/null || true
	-@cd $(OBJDIR) ; rm $(ACDC_D_OBJECTS) $(ACDC_D_MODS) 2>/dev/null || true
	-@cd $(OBJDIR) ; rm Superbox.o                       2>/dev/null || true
	-@rm superbox.exe                                    2>/dev/null || true

	# -@cd $(OBJDIR) ; rm $(UHMA_OBJECTS) $(UHMA_MODS)   2>/dev/null || true
# If you really want to remove chemistry objects too, use this
cleanall:
	-@cd $(OBJDIR) ; rm $(BOX_OBJECTS) $(BOX_MODS)       2>/dev/null || true
	-@cd $(OBJDIR) ; rm $(ACDC_OBJECTS) $(ACDC_MODS)      2>/dev/null || true
	-@cd $(OBJDIR) ; rm $(ACDC_D_OBJECTS) $(ACDC_D_MODS)  2>/dev/null || true
	-@cd $(OBJDIR) ; rm Superbox.o                        2>/dev/null || true
	-@rm superbox.exe                                     2>/dev/null || true
