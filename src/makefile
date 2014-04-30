# This makefile uses intel fortran compiler for linux

# Makefile variables -g 
FC      = ifort
FCFLAGS = -g  -parallel -mkl -traceback -CB
# flag -g increase runtime in raytracing pre-processing
# by a factor of almost 8
# FCFLAGS	 = -parallel -mkl -traceback -CB
VPATH	:= src
OBJDIR	:= ../obj
LIBPATH := /opt/intel/mkl/lib/intel64
LIBS    := $(LIBPATH)/libmkl_intel_lp64.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
  			$(LIBPATH)/libmkl_sequential.a \
  			$(LIBPATH)/libmkl_core.a  \
			$(LIBPATH)/libmkl_intel_ilp64.a \
  			-lpthread

# All compiled output
OBJECTS	:= $(OBJDIR)/test.o $(OBJDIR)/htfem.o \
		   $(OBJDIR)/readmesh.o $(OBJDIR)/element.o \
		   $(OBJDIR)/assembly.o $(OBJDIR)/boundary.o \
		   $(OBJDIR)/solver.o $(OBJDIR)/postproc.o \
		   $(OBJDIR)/transient.o
		   

# Make all target
all: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(OBJDIR)/test $(OBJECTS) $(LIBS)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/test.o: $(OBJDIR)/htfem.o $(OBJDIR)/readmesh.o \
				  $(OBJDIR)/element.o $(OBJDIR)/assembly.o \
				  $(OBJDIR)/boundary.o $(OBJDIR)/solver.o \
				  $(OBJDIR)/postproc.o $(OBJDIR)/transient.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c test.f90 -o $(OBJDIR)/test.o

$(OBJDIR)/htfem.o: $(OBJDIR)/readmesh.o $(OBJDIR)/element.o \
				   $(OBJDIR)/assembly.o $(OBJDIR)/boundary.o \
				   $(OBJDIR)/solver.o $(OBJDIR)/postproc.o \
				   $(OBJDIR)/transient.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c htfem.f90 -o $(OBJDIR)/htfem.o


$(OBJDIR)/transient.o: $(OBJDIR)/assembly.o $(OBJDIR)/solver.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c transient.f90 -o $(OBJDIR)/transient.o

$(OBJDIR)/postproc.o: $(OBJDIR)/element.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c postproc.f90 -o $(OBJDIR)/postproc.o

$(OBJDIR)/assembly.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c assembly.f90 -o $(OBJDIR)/assembly.o

$(OBJDIR)/boundary.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c boundary.f90 -o $(OBJDIR)/boundary.o

$(OBJDIR)/element.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c element.f90 -o $(OBJDIR)/element.o

$(OBJDIR)/readmesh.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c readmesh.f90 -o $(OBJDIR)/readmesh.o

$(OBJDIR)/solver.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c solver.f90 -o $(OBJDIR)/solver.o

# Raytracing objects
OBJECTS_rt := $(OBJDIR)/raytracing.o  $(OBJDIR)/rt_mods.o \
	      $(OBJDIR)/rt_math.o $(OBJDIR)/rt_preprocess.o \
	      $(OBJDIR)/rt_tracing.o $(OBJDIR)/rt_elements.o	   

# make only raytracing
rt: $(OBJECTS_rt) 
	$(FC) $(FCFLAGS) -o $(OBJDIR)/rt.out $(OBJECTS_rt) $(LIBS)

$(OBJDIR)/raytracing.o: $(OBJDIR)/rt_mods.o $(OBJDIR)/rt_math.o \
			$(OBJDIR)/rt_preprocess.o $(OBJDIR)/rt_tracing.o| $(OBJDIR) \
			$(OBJDIR)/rt_elements.o
	$(FC) $(FCFLAGS) -I$(OBJDIR) -c raytracing.f90 -o $(OBJDIR)/raytracing.o

$(OBJDIR)/rt_mods.o: | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt_mods.f90 -o $(OBJDIR)/rt_mods.o	
	
$(OBJDIR)/rt_math.o: $(OBJDIR)/rt_mods.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt_math.f90 -o $(OBJDIR)/rt_math.o

$(OBJDIR)/rt_elements.o: $(OBJDIR)/rt_math.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt_elements.f90 -o $(OBJDIR)/rt_elements.o
			
$(OBJDIR)/rt_preprocess.o: $(OBJDIR)/rt_elements.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt_preprocess.f90 -o $(OBJDIR)/rt_preprocess.o
	
$(OBJDIR)/rt_tracing.o: $(OBJDIR)/rt_elements.o | $(OBJDIR)
	$(FC) $(FCFLAGS) -I$(OBJDIR) -module $(OBJDIR) -c rt_tracing.f90 -o $(OBJDIR)/rt_tracing.o


# Cleaning everything
clean:
	rm -rf $(OBJDIR)
	rm -rf *.out
	rm -rf *.o *.mod test
