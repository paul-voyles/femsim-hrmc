HRMC = ../hrmc
FEMSIM = ../femsim

SRC = globals.f90 read_inputs.f90 model.f90 rmc_functions.f90 scattering_factors.f90 fem.f90 eam.f90 hrmc.f90
OBJ = $(addsuffix .o, $(basename $(SRC)))

LIBS=

FC = mpif90

.DEFAULT_GOAL := nothing

FCFLAGS = 

ifeq ($(MAKECMDGOALS),)
$(info )
$(info You must specify 'make hrmc' or 'make femsim')
$(info )
endif

# Set FCFLAGS for debugging if 'debug' is in targets, otherwise set normally
ifneq (,$(findstring debug,$(MAKECMDGOALS)))
FCFLAGS += -g -debug all -check all -implicitnone -warn all -openmp -fpp -DDEBUG
MAKECMDGOALS := $(filter-out debug,$(MAKECMDGOALS))
$(info Making debug version of $(MAKECMDGOALS)...)
else
FCFLAGS += -O2 -openmp -fpp
endif

ifneq (,$(findstring femsim,$(MAKECMDGOALS)))
FCFLAGS += -DFEMSIM
endif


.PHONY: clean debug nothing hrmc femsim
.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90
%.o: %.f90
	${FC} -c $< ${FCFLAGS}

hrmc: ${HRMC}

femsim: ${FEMSIM}

clean:
	rm -f *.mod *.o ../femsim ../hrmc

debug: ;
nothing: ;

${FEMSIM}: ${OBJ} hrmc.o
	${FC} ${OBJ} ${FCFLAGS} -o ${FEMSIM}
	rm hrmc.o  # Needed to remove hrmc.o so that it can be recompiled for hrmc correctly

${HRMC}: ${OBJ} hrmc.o
	${FC} ${OBJ} ${FCFLAGS} -o ${HRMC}
	rm hrmc.o  # Needed to remove hrmc.o so that it can be recompiled for femsim correctly
