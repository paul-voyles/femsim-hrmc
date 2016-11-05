MODULE HRMC_Global

  ! Define parameters used later in the calculation.
  double precision, parameter :: PI    = 3.141592653589793238462643383279502884
  double precision, parameter :: TWOPI = 6.283185307179586476925286766559005768
  integer, parameter :: ATOMS_PER_HUTCH = 1     ! average number of atoms in a hutch; best performance when this is 1
  double precision, parameter :: FEM_BIN_WIDTH = 0.005      ! bin width (in angstrom) of the gr_i in fem intensity calculation 
  integer, parameter ::  HRMC_STEPS = 1000      ! # of hrmc steps between model / sim data dumps
  ! MPI variables
  integer :: myid, numprocs, mpierr
  integer :: color, communicator, parent_comm
  character(len=16) :: color_str
  integer :: stdout = 6
  integer :: stderr = 0

END MODULE HRMC_Global
