! Hybrid Reverse Monte Carlo structural simulation

program hrmc

    use hrmc_global
    use readinputs
    use model_mod
    use fem_mod
    use hrmc_functions
    use eam_mod
    implicit none
    include 'mpif.h'
    ! HRMC / Femsim objects
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename
    character (len=256) :: eam_filename
    character (len=256) :: jobID, c, step_str
    character (len=256) :: vki_fn, vkf_fn, output_model_fn, final_model_fn, chi_squared_file, acceptance_rate_fn, femfile
    character (len=256) :: paramfile_restart
    real :: temperature
    real :: max_move
    real :: Q, res, alpha
    real, pointer, dimension(:) :: k
    double precision, pointer, dimension(:) :: vk, vk_exp, vk_exp_err, v_background
    real, allocatable, dimension(:,:) :: cutoff_r 
    real, pointer, dimension(:,:) :: scatfact_e
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: scale_fac, scale_fac_initial, beta, boltzmann
    double precision :: chi2_old, chi2_new, del_chi, chi2_gr, chi2_vk, chi2_no_energy, chi2_initial
    real :: R
    integer :: i, j, step_start, step_end
    integer :: w
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: istat, status2, length
    integer :: iseed2
    real :: randnum
    double precision :: te1, te2
    logical :: square_pixel, accepted
    integer :: ipvd, nthr
    real :: x! This is the parameter we will use to fit vsim to vas.
    integer, dimension(100) :: acceptance_array
    real :: avg_acceptance = 1.0
    integer :: temp_move_decrement
    character(3), dimension(118) :: syms

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

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

    if(myid .eq. 0) then
        call get_command_argument(1, c, length, istat)
        if (istat == 0) then
            jobID = "_"//trim(c)
        else
            error stop "istat for jobid get_command_arg was nonzero"
            !jobID = '_temp'
        end if
        call get_command_argument(2, c, length, istat)
        if (istat == 0) then
            param_filename = trim(c)
        else
            error stop "istat for paramfile get_command_arg was nonzero"
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
    call mpi_bcast(param_filename, 256, MPI_CHARACTER, 0, mpi_comm_world, mpierr)
    if(myid.eq.0) write(*,*) "Paramfile: ", trim(param_filename)
    
    ! Read input parameters
    call read_inputs(param_filename,model_filename, femfile, eam_filename, step_start, step_end, temp_move_decrement, temperature, max_move, cutoff_r, iseed2, alpha, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, scale_fac_initial, Q, status2)
    temperature = temperature*(sqrt(0.7)**(step_start/temp_move_decrement))
    max_move = max_move*(sqrt(0.94)**(step_start/temp_move_decrement))

    ! Read input model
    call read_model(model_filename, m, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    if(myid .eq. 0) then
    write(*,*) "Model filename: ", trim(model_filename)
    write(*,*)
    endif

    scale_fac = scale_fac_initial
    res = 0.61/Q
    nk = size(k)
    boltzmann = 8.6171e-05
    beta=1./((boltzmann)*temperature)
    square_pixel = .TRUE. ! HRMC uses square pixels, not round.

    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat,  square_pixel)
    allocate(vk(size(vk_exp)))
    vk = 0.0

    ! Print warning message if we are using too many cores.
    if(myid.eq.0) then
        call print_sampled_map(m, res, square_pixel)
        if(pa%npix /= 1) then
            if(numprocs > 3*nrot) write(0,*) "WARNING: You are using too many cores!"
        else
            if(numprocs > nrot) write(0,*) "WARNING: You are using too many cores!"
        endif
    endif

    !------------------- Call femsim. -----------------!

    ! Fem updates vk based on the intensity calculations and v_background.
    call fem(m, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)

    if(myid.eq.0)then
        ! Write initial vk 
        open(unit=52,file=trim(vki_fn),form='formatted',status='unknown')
            do i=1, nk
                write(52,*) k(i), vk(i)
            enddo
        close(52)
    endif

    !------------------- Call potential energy function. -----------------!

#ifndef FEMSIM
    call read_eam(m,eam_filename)
    call eam_initial(m,te1)
    te1 = te1/m%natoms
    if(myid .eq. 0) write(*,*) "Energy = ", te1

    !------------------- Start HRMC. -----------------!

    call mpi_barrier(mpi_comm_world, mpierr)

        ! Calculate initial chi2
        chi2_no_energy = chi_square(alpha, vk_exp, vk_exp_err, vk, scale_fac, nk)

        chi2_initial = chi2_no_energy
        chi2_old = chi2_no_energy + te1

        i=step_start
        if(myid.eq.0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step start =       ", i
            write(*,*) "   Step end =         ", step_end
            write(*,*) "   Decrement # =      ", temp_move_decrement
            write(*,*) "   Energy =           ", te1
            write(*,*) "   LSqF V(k) =        ", chi2_no_energy
            write(*,*) "   Temperature =      ", temperature
            write(*,*) "   Max Move=          ", max_move
            write(*,*)
            ! Reset energy/chi_squared file
            open(36,file=trim(chi_squared_file), form='formatted', status='unknown')
                write(36,*) "step, chi2, energy"
                write(36,*) i, chi2_no_energy, te1
            close(36)
            open(37,file=trim(acceptance_rate_fn), form='formatted', status='unknown', access='append')
                write(37,*) "step, acceptance rate averaged over last 1000 steps"
            close(37)
        endif


        ! HRMC loop begins. The loop never stops.
        do while (i .le. step_end)

            i = i+1
            if(myid .eq. 0) write(*,*) "Starting step", i

            call random_move(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, max_move)
            ! check_curoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m,cutoff_r,w) )
                ! Check_cutoffs returned false so reset positions and try again.
                m%xx%ind(w) = xx_cur
                m%yy%ind(w) = yy_cur
                m%zz%ind(w) = zz_cur
                call random_move(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, max_move)
            end do
            ! Update hutches, data for chi2, and chi2/del_chi
            call hutch_move_atom(m, w, xx_new, yy_new, zz_new)
    
            call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
            te2 = te2/m%natoms

            ! Calculate a randnum for accept/reject
            randnum = ran2(iseed2)

            call fem_update(m, w, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)

            chi2_no_energy = chi_square(alpha, vk_exp, vk_exp_err, vk, scale_fac, nk)
            chi2_new = chi2_no_energy + te2
            del_chi = chi2_new - chi2_old
            call mpi_bcast(del_chi, 1, mpi_double, 0, mpi_comm_world, mpierr)

            if(myid .eq. 0) write(*,*) "Energy = ", te2
            if(myid .eq. 0) write(*,*) "Del-V(k) = ", chi2_no_energy
            if(myid .eq. 0) write(*,*) "Del-chi = ", del_chi
            if(myid .eq. 0) write(*,*) "chi2_old = ", chi2_old
            if(myid .eq. 0) write(*,*) "chi2_new = ", chi2_new

            ! Test if the move should be accepted or rejected based on del_chi
            if(del_chi < 0.0) then
                ! Accept the move
                call fem_accept_move(mpi_comm_world)
                e1 = e2 ! eam
                chi2_old = chi2_new
                accepted = .true.
                if(myid.eq.0) write(*,*) "MC move accepted outright."
            else
                ! Based on the random number above, even if del_chi is negative, decide
                ! whether to move or not (statistically).
                if(log(1.-randnum) < -del_chi*beta) then
                    ! Accept move
                    call fem_accept_move(mpi_comm_world)
                    e1 = e2 ! eam
                    chi2_old = chi2_new
                    accepted = .true.
                    if(myid.eq.0) write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
                else
                    ! Reject move
                    call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                    call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                    call fem_reject_move(m, w, xx_cur, yy_cur, zz_cur, mpi_comm_world)
                    e2 = e1 ! eam
                    accepted = .false.
                    if(myid.eq.0) write(*,*) "MC move rejected."
                endif
            endif

            if(myid .eq. 0) then
                if(accepted) then
                    acceptance_array(mod(i,100)+1) = 1
                else
                    acceptance_array(mod(i,100)+1) = 0
                endif
                avg_acceptance = sum(acceptance_array)/100.0
                ! Writing to 0 is stderr
                !if(i .ge. 100 .and. avg_acceptance .le. 0.05 .and. mod(i,100) .eq. 0) write(0,*) "WARNING!  Acceptance rate is low:", avg_acceptance
            endif

            ! Periodically save data.
            if(myid .eq. 0) then
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
                if(accepted) then
                    ! Write chi2 and energy info
                    open(36,file=trim(chi_squared_file),form='formatted',status='unknown',access='append')
                        write(36,*) i, chi2_no_energy, te2
                    close(36)
                endif
                if(mod(i,100)==0 .and. i .ge. 100)then
                    ! Write to acceptance rate
                    open(40,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                        write(40,*) i, avg_acceptance
                    close(40)
                endif
            endif

            ! Every 'temp_move_decrement' steps lower the temp, max_move, and reset beta.
            if(mod(i,temp_move_decrement) .eq. 0 .and. i .ne. step_start)then
                temperature = temperature * sqrt(0.7)
                if(myid .eq. 0) write(*,*) "Lowering temp to", temperature, "at step", i
                max_move = max_move * sqrt(0.94)
                if(myid .eq. 0) write(*,*) "Lowering max_move to", max_move, "at step", i
                beta = 1./((boltzmann)*temperature)

                ! Write to param_resume.in file
                if(myid .eq. 0) then
                    open(unit=53,file=trim(paramfile_restart),form='formatted',status='unknown')
                        write(53,*) '# HRMC parameter file, generated to restart a sim ', trim(jobid(2:))
                        write(53,*) trim(output_model_fn)
                        write(53,*) trim(femfile)
                        write(53,*) trim(eam_filename)
                        write(53,*) step_start+i, step_end
                        write(53,*) temperature, max_move, temp_move_decrement
                        write(53,*) nelements
                        do j=1,nelements
                            write(53,*) cutoff_r(:,j)
                        enddo
                        write(53,*) iseed2
                        write(53,*) alpha
                        write(53,*) scale_fac_initial
                        write(53,*) nphi, npsi, ntheta
                        write(53,*) Q
                    close(53)
                endif
            endif

        enddo !HRMC do loop

        ! The hrmc loop finished. Write final data.
        if(myid .eq. 0) then
            write(*,*) "Monte Carlo Finished!"

            ! Write to param_resume.in file
            open(unit=53,file=trim(paramfile_restart),form='formatted',status='unknown')
                write(53,*) '# HRMC parameter file, generated to restart a sim ', trim(jobid(2:))
                write(53,*) trim(final_model_fn)
                write(53,*) trim(femfile)
                write(53,*) trim(eam_filename)
                write(53,*) step_start+i-1, step_end+i-1
                write(53,*) temperature, max_move, temp_move_decrement
                write(53,*) nelements
                do j=1,nelements
                    write(53,*) cutoff_r(:,j)
                enddo
                write(53,*) iseed2
                write(53,*) alpha
                write(53,*) scale_fac_initial
                write(53,*) nphi, npsi, ntheta
                write(53,*) Q
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

    call mpi_finalize(mpierr)

end program hrmc
