MODULE HRMC_Global

  ! Define parameters used later in the calculation.  These don't
  ! really need to be global, but it makes them easier to find.
  real, parameter :: PI = 3.141593
  real, parameter :: TWOPI = 6.283186
  integer, parameter :: ATOMS_PER_HUTCH = 1     ! average number of atoms in a hutch
  real, parameter :: FEM_BIN_WIDTH = 0.005      ! bin width (in angstrom) of the gr_i in fem intensity calculation 
  integer, parameter ::  HRMC_STEPS = 10000     ! # of hrmc steps between model / sim data dumps
  integer :: myid, numprocs, mpierr             ! mpi variables	

END MODULE HRMC_Global
