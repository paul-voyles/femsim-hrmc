!This module defines global variables
!Written by Feng Yi on 12/18/2008
! parameters added 12/20/08, pmv
! box size variables added 12/20/08 pmv
! add variable ATOMS_PER_HUTCH and an array variable atom_type, fy,12/23/2008
! changed composition to be 1D array to match atom_type, pmv, 12/29/08
! removed all model parameters.  they are now in derived type model 01-05-09 pmv

!*********************************************
!*****Define atomic information***************
!*********************************************

MODULE HRMC_Global

  ! Define parameters used later in the calculation.  These don't
  ! really need to be global, but it makes them easier to find.
  real, parameter :: PI = 3.141593
  real, parameter :: TWOPI = 6.283186
  integer, parameter :: ATOMS_PER_HUTCH = 1      ! average number of atoms in a hutch
  real, parameter :: FEM_BIN_WIDTH = 0.005       ! bin width (in angstrom) of the gr_i in fem intensity calculation 
  integer, parameter ::  HRMC_STEPS = 10000       ! # of hrmc steps between model / sim data dumps
                                                 ! calculation of simulated data
  integer :: myid, numprocs, mpierr       ! mpi variables	

END MODULE HRMC_Global
