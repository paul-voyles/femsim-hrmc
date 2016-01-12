APP = ../hrmc

SRC = ran2.f90 globals.f90 read_inputs.f90 model.f90 rmc_functions.f90 scattering_factors.f90 fem.f90 gr.f90 eam_fs.f90 hrmc.f90
OBJ = $(addsuffix .o, $(basename $(SRC)))

LIBS=

FC = mpif90
CC = icc


.PHONY: opt debug clean
.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90
%.o: %.f90
	${FC} -c $< ${FCFLAGS} ${LMP_INC}

# Opt is at the top to make it the default
# Need to define compiler flags for opt / debug
opt: FCFLAGS = -O2 -openmp -fpp $(FEMSIM)

debug: FCFLAGS = -g -debug all -check all -implicitnone -warn all -openmp -fpp $(FEMSIM)

${APP}: ${OBJ}
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

clean:
	rm -f *.mod *.o ${APP}
