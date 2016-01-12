
! femsim, completely written from scratch in Fortran, based
! on modules developed for hrmc.
!
! usage is:
!
!     femsim <model file>

program femsim

  use model_mod
  use fem_mod
  use readinputs

  implicit none
  include "mpif.h"

  type(model) :: m
  real res
  real, dimension(:,:), pointer :: scatfac_e
  character(len=256) :: c, model_filename, outbase, outfile
  integer :: len, istat, nk, i
  integer :: nphi, ntheta, npsi
  LOGICAL use_femsim
  ! Added by Jason from hrmc_v_1_0.f90 for read_input function
  real, pointer, dimension(:) :: vk_exp, vk_exp_err, v_background
  real :: temperature, max_move, Q, pixel_distance, scale_fac
  real, pointer, dimension(:,:) :: cutoff_r
  logical, dimension(4) :: used_data_sets
  real, dimension(4) :: weights
  real, pointer, dimension(:) :: gr_e, r_e, gr_e_err, gr_n, r_n, gr_x, r_x
  real :: rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x
  integer :: fem_algorithm, total_steps, status2
  logical :: square_pixel
  character (len=256):: param_filename
  real, pointer, dimension(:) :: k, vk
  real, allocatable, dimension(:) :: i_k

  write(*,*)
  write(*,*) "This is the dev version of Femsim!"
  write(*,*)

  ! Added by jason so I could pass mpi_comm_world to fem(...) so that I could
  ! pass the correct arguments. We want to be able to use the same fem function
  ! for HRMC and femsim, and HRMC passes in this argument so I made femsim do it
  ! too in order to hopefully get it to work correctly.
  call mpi_init(mpierr)
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

  ! parse command line for model name
  if (command_argument_count() /= 1) then
     write (*,*) 'femsim invoked with the wrong number of arguments.  Proper usage is:'
     write (*,*) ' '
     write (*,*) '    femsim <input model file>'
     write (*,*) ' '
     stop
  else
     call get_command_argument(1, c, len, istat)
     if (istat == 0) then
        model_filename = trim(c)
     else
        write (*,*) 'Cannot read model file name.  Exiting.'
        stop
     end if
  endif

  write (*,*) 'Fortran femsim v 1.0; 03-06-09'
  write (*,*) ' '

  call read_model(model_filename, c, m, istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open model file: ', model_filename
     write (*,*) 'Exiting.'
     stop
  endif
  ! Added by jason on 06/28/13
  call check_model(m, istat)
  ! recenter_model is called in read_model
  !call recenter_model(0.0, 0.0, 0.0, m)

  ! Added by jason to make this more similar to HRMC.
  param_filename = 'param_file.in'
  allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
  call read_inputs(param_filename,temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err,     gr_n, r_n, gr_x, r_x, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

  ! Added by jason so that I didn't have to input these. The new read_inputs
  ! reads in most of this.
  res = 0.61/Q
  nk = size(k)
  square_pixel = .FALSE. ! Femsim uses round pixels.
  use_femsim = .TRUE.
  outbase = "output/"//model_filename(1:size(model_filename)-4)//"_femsim"
  ! Set-up parameters and arrays for fem_initialize

  allocate(vk(size(vk_exp)), i_k(nk),stat=istat)
  if( istat /= 0) then
     write (*,*) 'Cannot allcoate memory for V(k) in top level program.'
     write (*,*) 'Exiting.'
     stop
  endif

  ! Open the V(k) output file
  outfile = trim(outbase) // '_vk.out'
  open(unit=300, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif

  ! Check to make sure we can open the intensities output file before doing
  ! the full calculation
  outfile = trim(outbase) // '_int.out'
  open(unit=301, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif
  close(301)

  ! Fem initialization
  write (*,*) ' '
  write (*,*) 'Initializing FEM calculation.'
  call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfac_e, istat, square_pixel)
  if (istat /= 0) then
     write (*,*) 'Failure in fem_initialize.'
     write (*,*) 'Exiting.'
     stop
  endif

  ! Fem calculate with hutch
  use_femsim = .TRUE.
  write (*,*) ' '
  write (*,*) 'Calculating V(k)'
  call fem(m, res, k, vk, v_background, scatfac_e, mpi_comm_world, istat, square_pixel, use_femsim)
  if (istat /=0 ) then
     write (*,*) 'Failure in fem subroutine.'
     write (*,*) 'Exiting.'
     stop
  endif

  ! Write the output files
  write (*,*) ' '
  write (*,*) 'Writing output to files.'
  
  write (300, *) 'k   V(k)'
  do i=1, nk
     write (*, *) k(i), vk(i)
     write (300, *) k(i), vk(i)
  enddo
  close(unit=300)

  call write_intensities(outfile, k, istat)
  if( istat /= 0) then
     write (*,*) 'Intensities file output failed to file: ',outfile
  endif

  
  close(301)
 

  call I_average(i_k)
  
  outfile = trim(outbase) // '_average_i.out'
  open(unit=302, file=outfile, form='formatted', status='replace', iostat=istat)
  if (istat /= 0) then
     write (*,*) 'Cannot open output file: ',outfile
     write (*,*) 'Exiting.'
     stop
  endif
  write(302,*) 'k i(k)'
  do i=1, nk
    write(*,*) k(i), i_k(i) !debug
    write(302,*) k(i), i_k(i)
  enddo

  close(302)


end program femsim
