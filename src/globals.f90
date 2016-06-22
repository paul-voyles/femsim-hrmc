MODULE HRMC_Global

  ! Define parameters used later in the calculation.  These don't
  ! really need to be global, but it makes them easier to find.
  real, parameter :: PI = 3.141593
  real, parameter :: TWOPI = 6.283186
  integer, parameter :: ATOMS_PER_HUTCH = 1     ! average number of atoms in a hutch; best performance when this is 1
  real, parameter :: FEM_BIN_WIDTH = 0.005      ! bin width (in angstrom) of the gr_i in fem intensity calculation 
  integer, parameter ::  HRMC_STEPS = 1000      ! # of hrmc steps between model / sim data dumps
  ! MPI variables
  integer :: myid, numprocs, mpierr
  integer :: color, communicator
  character(len=16) :: color_str

END MODULE HRMC_Global
