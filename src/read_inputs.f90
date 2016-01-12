!
!This module loads the experiment data.
!It was written by Jason Maldonis from the old read_inputs.f90 file on
!06/28/2013
!
module ReadInputs

contains

    subroutine read_inputs(param_filename, model_fn, femfile, eam_file, step_start, step_end, temp_move_decrement, temperature, max_move, cutoff_r, iseed2, alpha, V, k, V_err, V_background, ntheta, nphi, npsi, scale_fac, Q, status2)

    !param_filename=input file name containing initilizing parameters
    !temperature=beginning temperature for HRMC
    !max_move=maximum allowable movement distance
    !cutoff_r=cut off distance for different pairs
    !alpha = chi2 weighting factor
    !V=FEM intensity variance
    !k=scattering vector
    !V_err=FEM measurement variance error
    !V_background=
    !ntheta, nphi, npsi=rotation number for calculating simulated V
    !scale_fac= 1/3ts/te?
    !Q=resolution for FEM
    !fem_algorithm=which algorithm is used to do variance calculation
    !status2=whether param_filename is opened with success or not

        implicit none
        character (len=256), intent(in) :: param_filename  !Assume size array, be careful
        character (len=256), intent(out) :: model_fn, eam_file
        character (len=256), intent(inout) :: femfile ! The fem data file name, at most 20 characters
        integer, intent(out) :: step_start, step_end, temp_move_decrement
        real, intent(out) :: temperature
        real, intent(out) :: max_move
        real, intent(out), dimension(:,:), allocatable :: cutoff_r
        integer, intent(out) :: iseed2
        real, intent(out) :: alpha
        real, pointer, dimension(:) :: k
        double precision, pointer, dimension(:) :: v, v_err, v_background
        integer, intent(out) :: ntheta, nphi, npsi
        real, intent(out) :: q
        real, intent(out) :: scale_fac
        integer, intent(out) :: status2
        integer :: nelements

        !Variables declared in this subroutine, local variables
        !comment1=comment in param_filename
        character (len=256) comment1  
        character (len=256) scatteringfile ! The scattering file name, at most 20 characters
        integer filenamelength !The file name length in scatteringfile or femfile
        integer i
        real indicator_end !Indicator_end=-1 means the end of file reaches
        integer status1 !Indicate the status of opening file in this subroutine
        logical file_end !Indicate fileend has reached
            ! Files for electron, neutron, x-ray scattering and fem data
        integer num_line  !Number of lines in each data file except comment line
        real, pointer, dimension(:) :: tempdata !Temperature data
        integer stat_allocate1, stat_allocate2, stat_allocate3,stat_allocate4 !Allocate status, 0 means success

        open(20, file=param_filename,iostat=status2, status='old')
        if(status2 .ne. 0) then !open fails
            print *, 'cannot open file with name: ', param_filename
            return
        endif
        read(20, '(A256)') comment1 !read comment and it is neglected in the input
        read(20, '(A256)') model_fn; model_fn = adjustl(model_fn)
        read(20, '(A256)') femfile; femfile= adjustl(femfile)
        read(20, '(A256)') eam_file; eam_file= adjustl(eam_file)
        read(20, *) step_start, step_end
        read(20, *) temperature, max_move, temp_move_decrement
        read(20, * ) nelements
        allocate(cutoff_r(nelements,nelements))
        read(20, *) cutoff_r ! This is a nelements by nelements matrix
        read(20, *) iseed2
        read(20, *) alpha
        read(20, *) scale_fac
        read(20, *) nphi, npsi, ntheta
        read(20, *) q
        close(20)

        num_line = 0
        filenamelength=len_trim(femfile)
        file_end=.false.
        open(30,file=femfile(1:filenamelength),iostat=status1,status='old') 
        if(status1 .eq. 0) then !open succeeds
            read(30, '(a80)') comment1 !first line is comment
            ! Count how many data pairs are in the file. -1 denotes EOF
            do while( .not. file_end)
                read(30, *) indicator_end
                if(abs(indicator_end+1.0) .le. 1e-6) then !file end is reached
                    exit !exit this loop
                else
                    num_line = num_line + 1
                endif
            enddo
            rewind(30) !go to the beginning of the file
            read(30, '(a80)') comment1
            allocate(tempdata(4*num_line),stat=stat_allocate1)
            !read k ,v, and v_err data
            !read k first, then v, last v_err
            if(stat_allocate1 .eq. 0) then
                read(30, *) tempdata
                allocate(k(num_line), stat=stat_allocate2)
                allocate(v(num_line), stat=stat_allocate3)
                allocate(v_err(num_line), stat=stat_allocate4)
                allocate(v_background(num_line))

                if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0) .and. (stat_allocate4 .eq. 0)) then
                    k=tempdata(1:4*num_line:4)
                    v=tempdata(2:4*num_line:4)
                    v_err=tempdata(3:4*num_line:4)
                    v_background=tempdata(4:4*num_line:4)
                else
                    print *, 'fem part 2 or 3, or 4 fails!'
                    return
                endif !allocate2, 3 and 4
            else
                print *, 'fem part allocation fail'
                return
            endif
        deallocate(tempdata)
        else
            print *, 'open fem file fails', femfile(1:filenamelength)
        endif
        close(30)
    end subroutine read_inputs

end module readinputs

