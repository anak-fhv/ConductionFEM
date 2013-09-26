# This makefile uses intel fortran compiler for linux

# Makefile variables
FC		:= ifort
FCFLAGS	:= -g
VPATH	:= src
OBJDIR	:= ../obj

# All compiled output
OBJECTS	:= $(OBJDIR)/test.o $(OBJDIR)/finitesolver.o $(OBJDIR)/readmesh.o $(OBJDIR)/numericalmethods.o

# Make all target
all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(OBJDIR)/test $(OBJECTS)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/test.o: $(OBJDIR)/finitesolver.o $(OBJDIR)/readmesh.o $(OBJDIR)/numericalmethods.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c test.f90 -o $(OBJDIR)/test.o

$(OBJDIR)/finitesolver.o: $(OBJDIR)/readmesh.o $(OBJDIR)/numericalmethods.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c finitesolver.f90 -o $(OBJDIR)/finitesolver.o

$(OBJDIR)/readmesh.o: $(OBJDIR)/numericalmethods.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c readmesh.f90 -o $(OBJDIR)/readmesh.o

$(OBJDIR)/numericalmethods.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c numericalmethods.f90 -o $(OBJDIR)/numericalmethods.o

# Cleaning everything
clean:	
	rm -rf $(OBJDIR)
	rm -rf *.o *.mod test

# End of the makefile
