# This makefile uses intel fortran compiler for linux

# Makefile variables
FC		 = ifort
FCFLAGS	 = -g -openmp -mkl -traceback
VPATH	:= src
OBJDIR	:= ../obj
LIBPATH := /opt/intel/mkl/lib/intel64
LIBS := 	$(LIBPATH)/libmkl_intel_lp64.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
			$(LIBPATH)/libmkl_intel_ilp64.a \
  			-lpthread

# All compiled output
OBJECTS	:= $(OBJDIR)/test.o $(OBJDIR)/conductionfem.o $(OBJDIR)/readmesh.o $(OBJDIR)/elementcalculations.o

# Make all target
all: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(OBJDIR)/test $(OBJECTS) $(LIBS)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/test.o: $(OBJDIR)/conductionfem.o $(OBJDIR)/readmesh.o $(OBJDIR)/elementcalculations.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c test.f90 -o $(OBJDIR)/test.o 

$(OBJDIR)/conductionfem.o: $(OBJDIR)/readmesh.o $(OBJDIR)/elementcalculations.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c conductionfem.f90 -o $(OBJDIR)/conductionfem.o

$(OBJDIR)/readmesh.o: $(OBJDIR)/elementcalculations.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c readmesh.f90 -o $(OBJDIR)/readmesh.o 

$(OBJDIR)/elementcalculations.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c elementcalculations.f90 -o $(OBJDIR)/elementcalculations.o

# Cleaning everything
clean:	
	rm -rf $(OBJDIR)
	rm -rf *.out
	rm -rf *.o *.mod test

# End of the makefile
