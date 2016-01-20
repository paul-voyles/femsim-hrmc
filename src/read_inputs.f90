!This module loads the experiment data.

module ReadInputs

contains

    subroutine read_inputs(param_filename, femsim, model_fn, femfile, Q, ntheta, nphi, npsi, scale_fac, eam_file, step_start, step_end, temp_move_decrement, temperature, max_move, cutoff_r, iseed2, alpha, V, k, V_err)
        implicit none
        ! Inputs
        character (len=256), intent(in) :: param_filename  ! Assume size array, be careful
        logical, intent(in) :: femsim
        character (len=256), intent(out) :: model_fn, eam_file, femfile
        integer, intent(out) :: step_start, step_end, temp_move_decrement
        real, intent(out) :: temperature
        real, intent(out) :: max_move
        real, intent(out), dimension(:,:), allocatable :: cutoff_r
        integer, intent(out) :: iseed2
        real, intent(out) :: alpha
        real, pointer, dimension(:) :: k
        double precision, pointer, dimension(:) :: v, v_err
        integer, intent(out) :: ntheta, nphi, npsi
        real, intent(out) :: q
        real, intent(out) :: scale_fac

        ! Local variables
        character (len=256) :: comment
        character (len=256) :: line
        integer :: status2
        integer :: nelements
        integer :: filenamelength ! The file name length in scatteringfile or femfile
        integer :: i, line_size ! Counter i; number of numbers on a line
        integer :: status1 ! Indicate the status of opening file in this subroutine
        integer :: num_line  ! Number of lines in each data file except comment line

        open(20, file=param_filename,iostat=status2, status='old')
        if(status2 .ne. 0) then !open fails
            print *, 'cannot open file with name: ', param_filename
            return
        endif
        read(20, '(A256)') comment ! Read comment; it is ignored
        read(20, '(A256)') model_fn; model_fn = adjustl(model_fn)
        read(20, '(A256)') femfile; femfile = adjustl(femfile)
        read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
        read(line, *) q
        read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
        read(line, *) nphi, npsi, ntheta
        read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
        read(line, *) scale_fac
        if(.not. femsim) then
            read(20, '(A256)') eam_file; eam_file= adjustl(eam_file)
            read(20, '(A)') line
            read(line, *) step_start, step_end
            read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
            read(line, *) temperature, max_move, temp_move_decrement
            read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
            read(line, *) iseed2
            read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
            read(line, *) alpha
            read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
            read(line, * ) nelements
            allocate(cutoff_r(nelements,nelements))
            do i=1, nelements
                read(20, '(A)') line; if(index(line, "#") .ne. 0) then; line = line(1:index(line, "#")-1); line = adjustl(line); endif
                read(line, *) cutoff_r(:,i) ! This is a nelements by nelements matrix
            end do
        endif
        close(20)

        if(index(model_fn, "#") .ne. 0) then
            model_fn = model_fn(1:index(model_fn, "#")-1)
            model_fn = adjustl(model_fn)
        endif
        if(index(femfile, "#") .ne. 0) then
            femfile = femfile(1:index(femfile, "#")-1)
            femfile = adjustl(femfile)
        endif
        if(.not. femsim) then
            if(index(eam_file, "#") .ne. 0) then
                eam_file = eam_file(1:index(eam_file, "#")-1)
                eam_file = adjustl(eam_file)
            endif
        endif

        filenamelength = len_trim(femfile)
        open(30, file=femfile(1:filenamelength), iostat=status1, status='old') 
        if(status1 .eq. 0) then ! Open succeeds
            read(30, '(a80)') comment ! First line is comment
            num_line = 0
            do while(.true.) ! Until end of file
                read(30, '(A)', iostat=status1) line
                if(status1 .ne. 0) exit
                if(num_line .gt. 0) then
                    if(line_size .ne. count_numbers_on_line(line)) then
                        stop "Each line must have the same number of elements on it in your FEM data file!"
                    endif
                endif
                line_size = count_numbers_on_line(line)
                line = adjustl(line)
                num_line = num_line + 1
            enddo

            if(line_size .eq. 1 .and. .not. femsim) stop "V(k) data must be included in the FEM data file."
            if(line_size .ne. 1 .and. line_size .ne. 3) then
                write(0,*) "The number of columns in the FEM datafile must be 1 or 3 (got", line_size,")"
                stop
            endif

            if(line_size .ge. 1) then  ! k data only (for FEMSIM)
                allocate(k(num_line))
            endif
            if(line_size .ge. 3) then  ! k, vk, and vkerr
                allocate(v(num_line))
                allocate(v_err(num_line))
            endif

            rewind(30) ! Go to the beginning of the file
            read(30, '(a80)') comment

            do i=1, num_line
                if(line_size .eq. 1) then  ! k data only (for FEMSIM)
                    read(30, *) k(i)
                else if(line_size .eq. 3) then  ! k, vk, and vkerr
                    read(30, *) k(i), v(i), v_err(i)
                endif
            end do
        else
            write(0,*) "Failed to open the FEM data file ", femfile(1:filenamelength)
        endif
        close(30)
    end subroutine read_inputs

    function count_numbers_on_line(str) result(nums)
        character(*) :: str
        logical :: sep
        integer :: i, nums
        sep = .true.
        nums = 0
        do i = 1, len(str)
            if ( (ichar(str(i:i)) .ge. 48 .and. ichar(str(i:i)) .le. 57) .or.  ichar(str(i:i)) .eq. 46 ) then ! we have a number or a .
                if(sep) nums = nums + 1
                sep = .false.
            else
                sep = .true.
            end if
        end do
    end function

end module readinputs
