! Hybrid Reverse Monte Carlo using FEMSIM


program hrmc

    use hrmc_global
    use readinputs
    use model_mod
    use fem_mod
    use hrmc_functions
    use eam_mod
    implicit none
#ifndef SERIAL
    include 'mpif.h'
#else
    integer :: mpi_comm_world
#endif
    ! Femsim  variables
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename
    character (len=256) :: eam_filename
    character (len=256) :: jobID, c
    character (len=256) :: vki_fn, vkf_fn, final_model_fn, chi_squared_file, acceptance_rate_fn, femfile
    character (len=256) :: paramfile_restart
    real :: temperature
    real :: max_move
    double precision :: Q, res, alpha
    double precision, pointer, dimension(:) :: k, vk, vk_exp, vk_exp_err
    double precision, allocatable, dimension(:,:) :: cutoff_r 
    double precision, pointer, dimension(:,:) :: scatfact_e
    double precision :: scale_fac, scale_fac_initial, beta, boltzmann
    integer :: i
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: step_start, step_end
    integer :: istat, length
    integer :: seed
    logical :: femsim
    integer :: ipvd
    integer :: temp_move_decrement
    character(3), dimension(118) :: syms
#ifndef FEMSIM
    ! HRMC variables
    character (len=256) :: output_model_fn
    character (len=256) :: step_str
    double precision :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    double precision :: chi2_prev_step, chi2_old, chi2_new, del_chi, chi2_no_energy, chi2_initial
    integer :: atom, j
    double precision :: randnum
    double precision :: te1, te2
    logical :: accepted, energy_accepted
    integer, dimension(100) :: acceptance_array
    real :: avg_acceptance = 1.0
#endif

    syms = (/ "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",  &
        "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",    &
        "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",&
        "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", &
        "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", &
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",&
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", &
        "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", &
        "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",&
        "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo" /)

    !------------------- Program setup. -----------------!
#ifdef FEMSIM
    femsim = .true.
#else
    femsim = .false.
#endif

#ifndef SERIAL
    ! Read in this special environemnt variable set by OpenMPI
    ! Context: http://stackoverflow.com/questions/35924226/openmpi-mpmd-get-communication-size
    ! (See the solution by Hristo Iliev)
    ! The below two lines set `color` to be 0, 1, ..., P where P is the number
    ! of executables given on the command line within `mpiexec` that are colon-deliminated
    call get_environment_variable("OMPI_MCA_orte_app_num", color_str)
    call str2int(color_str, color, istat)

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr)
    ! Split the world communicator into separate pieces for each mpiexec subprogram / spawn multiple
    call mpi_barrier(mpi_comm_world, mpierr)
    call mpi_comm_split(mpi_comm_world, color, myid, communicator, mpierr)

    call mpi_comm_get_parent(parent_comm, mpierr)

    ! Get each core's rank within its new communicator
    call mpi_comm_rank(communicator, myid, mpierr)
    call mpi_comm_size(communicator, numprocs, mpierr)
    call get_environment_variable("OMPI_MCA_orte_app_num", color_str)
    if(myid .eq. 0) write(*,*) "Successfully initialized MPI_COMM_WORLD"
    write(*,'(A10, I2, A4, I2, A11, I2, A6, I2, A18, I2)') "I am core ", myid, " of ", numprocs, " with color", color, ", root", 0, ", and communicator", communicator


    call mpi_barrier(communicator, mpierr)
#else
    myid = 0
    communicator = 0
    numprocs = 1
#endif

    if(myid .eq. 0) then
        call get_command_argument(1, c, length, istat)
        if (istat == 0) then
            jobID = "_"//trim(c)
        else
            write(0, *) "istat for jobid get_command_arg was nonzero"
            stop 1
        end if
        call get_command_argument(2, c, length, istat)
        if (istat == 0) then
            param_filename = trim(c)
        else
            write(0, *) "istat for paramfile get_command_arg was nonzero"
            stop 1
        end if
        param_filename = trim(param_filename)

        ! Set output filenames.
        write(vki_fn, "(A10)") "vk_initial"
        vki_fn = trim(trim(vki_fn)//jobID)//".txt"
        write(vkf_fn, "(A8)") "vk_final"
        vkf_fn = trim(trim(vkf_fn)//jobID)//".txt"
        write(final_model_fn, "(A11)") "model_final"
        final_model_fn = trim(trim(final_model_fn)//jobID)//".xyz"
        write(chi_squared_file, "(A11)") "chi_squared"
        chi_squared_file = trim(trim(chi_squared_file)//jobID)//".txt"
        write(acceptance_rate_fn, "(A15)") "acceptance_rate"
        acceptance_rate_fn = trim(trim(acceptance_rate_fn)//jobID)//".txt"
        write(paramfile_restart, "(A16)") "param_resume.in"
    endif

    !------------------- Read inputs and initialize. -----------------!

    ! Set input filenames.
    if(myid.eq.0) write(*,*) "Paramfile: ", trim(param_filename)
#ifndef SERIAL
    call mpi_bcast(param_filename, 256, MPI_CHARACTER, 0, communicator, mpierr)
#endif
    
    ! Read input parameters
    call read_inputs(param_filename, femsim, model_filename, femfile, Q, ntheta, nphi, npsi, scale_fac_initial, eam_filename, step_start, step_end, temp_move_decrement, temperature, max_move, cutoff_r, seed, alpha, vk_exp, k, vk_exp_err)

    call random_seed(put=(/seed/))

    ! Read input model
    call read_model(model_filename, m, istat)
    call check_model(m, istat)
    call recenter_model(dble(0.0), dble(0.0), dble(0.0), m)

    if(myid .eq. 0) then
    write(*,*) "Model filename: ", trim(model_filename)
    write(*,*)
    endif

    scale_fac = scale_fac_initial
    res = 0.61/Q
    nk = size(k)

    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat)
    allocate(vk(size(vk_exp)))
    vk = 0.0

    ! Print warning message if the user is using too many cores.
    ! For example, in a 1 pixel model, the number of cores should
    ! never be larger than the number of rotations.
    if(myid.eq.0) then
        call print_sampled_map(m, res)
        if(pa%npix /= 1) then
            if(numprocs > 3*nrot) write(0,*) "WARNING: You are using too many cores!"
        else
            if(numprocs > nrot) write(0,*) "WARNING: You are using too many cores!"
        endif
    endif

    !------------------- Call femsim. -----------------!

    ! Fem updates vk based on the intensity calculations.
    call fem(m, res, k, vk, scatfact_e, communicator, istat)

    if(myid.eq.0)then
        ! Write initial vk to file
        open(unit=52,file=trim(vki_fn),form='formatted',status='unknown')
            do i=1, nk
                write(52,*) k(i), vk(i)
            enddo
        close(52)
    endif

#ifndef FEMSIM
    temperature = temperature*(sqrt(0.7)**(step_start/temp_move_decrement))
    max_move = max_move*(sqrt(0.94)**(step_start/temp_move_decrement))
    boltzmann = 8.6171e-05
    beta=1./((boltzmann)*temperature)

    !------------------- Call potential energy function. -----------------!

    call read_eam(m,eam_filename)
    call eam_initial(m,te1)
    te1 = te1/m%natoms

    !------------------- Start HRMC. -----------------!

#ifndef SERIAL
    call mpi_barrier(communicator, mpierr)
#endif

        ! Calculate initial chi2
        chi2_no_energy = chi_square(alpha, vk_exp, vk_exp_err, vk, scale_fac, nk)
        chi2_initial = chi2_no_energy
        chi2_old = chi2_no_energy + te1

        i = step_start
        if(myid .eq. 0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step start =       ", i
            write(*,*) "   Step end =         ", step_end
            write(*,*) "   Decrement # =      ", temp_move_decrement
            write(*,*) "   Energy =           ", te1
            write(*,*) "   LSqF V(k) =        ", chi2_no_energy
            write(*,*) "   Cost Function =    ", chi2_old
            write(*,*) "   Temperature =      ", temperature
            write(*,*) "   Max Move =         ", max_move
            write(*,*)
            ! Reset energy/chi_squared and acceptance rate files
            open(36,file=trim(chi_squared_file), form='formatted', status='unknown')
                write(36,*) "step, chi2, energy"
                write(36,*) i, chi2_no_energy, te1
            open(40,file=trim(acceptance_rate_fn), form='formatted', status='unknown')
                write(40,*) "step, acceptance rate averaged over last 1000 steps"
        endif


        ! HRMC loop begins.
        te2 = te1
        e2 = e1
        do while (i .le. step_end)
            i = i + 1
            if(myid .eq. 0) write(*,*) "Starting step", i

            chi2_prev_step = chi2_old
            call random_move(m, atom, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, max_move)
            ! check_cutoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m, cutoff_r, atom) )
                ! check_cutoffs returned false so reset positions and try again.
                m%xx%ind(atom) = xx_cur
                m%yy%ind(atom) = yy_cur
                m%zz%ind(atom) = zz_cur
                call random_move(m, atom, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, max_move)
            end do
            ! Update hutches with the moved atom position
            call hutch_move_atom(m, atom, xx_new, yy_new, zz_new)
    
            call eam_mc(m, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
            te2 = te2/m%natoms

            ! Calculate a randnum for accept/reject
            call random_number(randnum)
            randnum = 1.0-randnum

            ! Decide whether to reject just based on the energy.
            ! The fit to V(k) cannot be negative, so there is a lower
            ! bound on chi2_no_energy (zero). So set it to 0 and check
            ! if the energy is so bad that the no value of chi2_no_energy
            ! can make the change in the objective function good enough
            ! to be accepted.
            ! This can save an entire call to fem_update.
            energy_accepted = .true.
            if(log(randnum) > -(te2-chi2_old)*beta) then
                energy_accepted = .false.
                accepted = .false.
                call reject_position(m, atom, xx_cur, yy_cur, zz_cur)
                call hutch_move_atom(m, atom, xx_cur, yy_cur, zz_cur)
                e2 = e1
            endif

            ! If the energy is good enough, run FEMSIM and re-evaluate the acceptance
            if(energy_accepted) then
                call fem_update(m, atom, res, k, vk, scatfact_e, communicator, istat)

                chi2_no_energy = chi_square(alpha, vk_exp, vk_exp_err, vk, scale_fac, nk)
                chi2_new = chi2_no_energy + te2
                del_chi = chi2_new - chi2_old
#ifndef SERIAL
                call mpi_bcast(del_chi, 1, mpi_double, 0, communicator, mpierr)
#endif

                ! Test if the move should be accepted or rejected based on del_chi
                if(del_chi < 0.0) then
                    ! Accept move
                    call fem_accept_move()
                    e1 = e2
                    chi2_old = chi2_new
                    accepted = .true.
                else
                    ! Based on the random number above, even if del_chi is negative, decide
                    ! whether to accept or not.
                    if(log(randnum) < -del_chi*beta) then
                        ! Accept move
                        call fem_accept_move()
                        e1 = e2 ! eam
                        chi2_old = chi2_new
                        accepted = .true.
                    else
                        ! Reject move
                        call reject_position(m, atom, xx_cur, yy_cur, zz_cur)
                        call hutch_move_atom(m,atom,xx_cur, yy_cur, zz_cur)
                        call fem_reject_move(m, atom, xx_cur, yy_cur, zz_cur)
                        e2 = e1
                        accepted = .false.
                    endif
                endif
            endif

            ! With the exception of decrement the temperature and max_move,
            ! everything after this line is data output to save information to disk.

            if(myid .eq. 0) then
                if(accepted) then
                    ! Write chi2 and energy info
                    write(36,*) i, chi2_no_energy, te2

                    ! Update acceptance rate
                    acceptance_array(mod(i,100)+1) = 1

                    if(del_chi < 0.0) then
                        write(*,*) "MC move accepted outright."
                    else
                        write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
                    endif

                    write(*,*) "Energy = ", te2
                    write(*,*) "Del-V(k) = ", chi2_no_energy
                    write(*,*) "Del-chi = ", del_chi
                    write(*,*) "chi2_prev_step = ", chi2_prev_step
                    write(*,*) "chi2_new = ", chi2_new

                else if(.not. energy_accepted) then
                    write(*,*) "MC move rejected solely due to energy.", te2
                else
                    write(*,*) "MC move rejected.", chi2_old, chi2_new
                    acceptance_array(mod(i,100)+1) = 0
                endif

                ! Write to acceptance rate
                if(mod(i,100)==0 .and. i .ge. 100)then
                    avg_acceptance = sum(acceptance_array)/100.0
                    write(40,*) i, avg_acceptance
                endif

                ! Periodically save data.
                if(mod(i,HRMC_STEPS) .eq. 0 .and. i .gt. step_start) then
                    write(output_model_fn, "(A12)") "model_update"
                    write(step_str,*) i
                    output_model_fn = trim(trim(trim(trim(output_model_fn)//jobID)//"_")//adjustl(step_str))//".xyz"
                    open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                        write(33,*) m%natoms
                        write(33,'(3F20.14)')m%lx,m%ly,m%lz
                        do j=1,m%natoms
                            write(33,'(A4, 3F20.14)') syms(m%znum%ind(j)), m%xx%ind(j), m%yy%ind(j), m%zz%ind(j)
                        enddo
                    close(33)
                endif
            endif

            ! Every 'temp_move_decrement' steps lower the temp, max_move, and reset beta.
            if(mod(i,temp_move_decrement) .eq. 0 .and. i .ne. step_start)then
                temperature = temperature * sqrt(0.7)
                max_move = max_move * sqrt(0.94)
                beta = 1./((boltzmann)*temperature)

                ! Write to param_resume.in file
                if(myid .eq. 0) then
                    write(*,*) "Lowering temp to", temperature, "at step", i
                    write(*,*) "Lowering max_move to", max_move, "at step", i

                    open(unit=53,file=trim(paramfile_restart),form='formatted',status='unknown')
                        write(53,*) '# HRMC parameter file, generated to restart a sim ', trim(jobid(2:))
                        write(53,*) trim(output_model_fn)
                        write(53,*) trim(femfile)
                        write(53,*) Q
                        write(53,*) nphi, npsi, ntheta
                        write(53,*) scale_fac_initial
                        write(53,*) trim(eam_filename)
                        write(53,*) step_start+i, step_end
                        write(53,*) temperature, max_move, temp_move_decrement
                        write(53,*) seed
                        write(53,*) alpha
                        write(53,*) nelements
                        do j=1,nelements
                            write(53,*) cutoff_r(:,j)
                        enddo
                    close(53)
                endif
            endif

        enddo ! HRMC loop

        ! The HRMC loop finished. Write final data.
        if(myid .eq. 0) then
            write(*,*) "HRMC Finished!"

            ! Close the energy/chi_squared and acceptance rate files
            close(36)
            close(40)

            ! Write to param_resume.in file
            open(unit=53,file=trim(paramfile_restart),form='formatted',status='unknown')
                write(53,*) '# HRMC parameter file, generated to restart a sim ', trim(jobid(2:))
                write(53,*) trim(final_model_fn)
                write(53,*) trim(femfile)
                write(53,*) Q
                write(53,*) nphi, npsi, ntheta
                write(53,*) scale_fac_initial
                write(53,*) trim(eam_filename)
                write(53,*) step_start+i-1, step_end+i-1
                write(53,*) temperature, max_move, temp_move_decrement
                write(53,*) seed
                write(53,*) alpha
                write(53,*) nelements
                do j=1,nelements
                    write(53,*) cutoff_r(:,j)
                enddo
            close(53)

            ! Write final vk
            open(unit=54,file=trim(vkf_fn),form='formatted',status='unknown')
            do i=1, nk
                write(54,*)k(i),vk(i)
            enddo
            close(54)

            ! Write final model
            open(unit=55,file=trim(final_model_fn),form='formatted',status='unknown')
            write(55,*) m%natoms
            write(55,'(3F20.14)')m%lx,m%ly,m%lz
            do i=1,m%natoms
                write(55,'(A4, 3F20.14)') syms(m%znum%ind(i)), m%xx%ind(i), m%yy%ind(i), m%zz%ind(i)
            enddo
        endif
#endif

    ! Free the sub-communicators and finalize
    if(parent_comm .ne. mpi_comm_null) then
        call mpi_comm_disconnect(parent_comm, mpierr)
        write(*,*) "Disconnected from parent!"
    endif
    call mpi_finalize(mpierr)

    write(*,*) "Successfully finished FEMSIM"
end program hrmc


elemental subroutine str2int(str, int, stat)
    implicit none
    character(len=*), intent(in) :: str
    integer, intent(out)         :: int
    integer, intent(out)         :: stat

    read(str, *, iostat=stat)  int
end subroutine str2int
