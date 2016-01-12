!*********************************************************
!**********This module is written by FY on 12/19/2008****
!**********Include three functions***********************
!*********************************************************

!    function chi_square, check_cutoffs and generate_move are included
!


MODULE hrmc_functions

  use model_mod
  implicit none

  interface
     function ran2(idum)
       real :: ran2
       integer :: idum
     end function ran2
  end interface

CONTAINS

    !***************************************************************
    !This function is used to calculate chi square
    function chi_square(alpha,vk_exp, vk_exp_err, vk, scale_fac, nk)

        !V=FEM intensity variance from input file
        !V_err=FEM measurement variance error
        !V_sim=simulated V(k) for variance data
        real :: alpha
        double precision, pointer, dimension(:) :: vk_exp, vk_exp_err
        double precision, pointer, dimension(:) :: vk
        real, intent(in) :: scale_fac
        integer, intent(in) ::nk
        integer i, j
        integer nf   !normalization factor   - jwh 04/25/2009
        double precision :: chi_square
        chi_square=0.0
        
        !FEM data
        nf = 0
        do j=1,nk
            nf = nf + 1
            !write(*,*) "DEBUG: v(j)*scale_fac-v_sim(j)=", v(j)*scale_fac-v_sim(j)
            !write(*,*) "DEBUG: scale_fac*V_err(j)=", scale_fac*V_err(j)
            !write(*,*) "DEBUG: Total addition is:", alpha*((v(j)*scale_fac-v_sim(j))/(scale_fac*V_err(j)))**2
            chi_square = chi_square+alpha*((vk_exp(j)*scale_fac-vk(j))/(scale_fac*vk_exp_err(j)))**2
        enddo
        chi_square = chi_square/nf

    end function chi_square


    subroutine update_scale_factor(sf, sf_initial, vsim, vas)
    ! Compute chi2(vsim*x, vas) for variable x to minimize chi2.
    ! (Note, I mean "variable" as in "changing", not the noun.)
    ! Set sf = sf/x.
    ! This isn't going to work because calling this multiple times will
    ! eventually monotonically change sf because vsim vs. vas will always be a
    ! little smaller or a little bigger (likely) and therefore sf will get worse
    ! and worse after teh first call. I need another variable I think.
        real, intent(out) :: sf
        real, intent(in) :: sf_initial
        ! vsim is the simulated vk as computed by the fast intensity algorithm. 
        ! vas is the simulated vk as computed by autoslice.
        double precision, intent(in), pointer, dimension(:) :: vsim, vas 
        real :: x! This is the parameter we will use to fit vsim to vas.
        integer :: i

        ! Below is the numerical solved solution to mimize
        ! chi2 = summation( (vsim(i)*x - vas(i))**2/vas(i), i=1, i=size(vsim) )
        ! by varying x. x will then be the weighting factor that gets vsim as
        ! close to vas as possible. I.e. vsim*x ~= vas.
        x = 0
        do i=1, size(vsim)
            x = x + vas(i)/vsim(i)
        enddo

        sf = sf_initial/x ! We need to divide because vk_exp ~= vsim/sf
        write(*,*) "Scale factor updated. Old,New =", sf_initial, sf
        write(*,*) "Conversion factor x =", x
    end subroutine update_scale_factor



!SUBROUTINE Calc_Type_Order(m, Type_Order) must be called first, then 
!this function can be called
!The check_cufoffs function return .TRUE. if 
!****no two atoms are within cutoff distance
!Otherwise return .FALSE.

    function check_cutoffs(m,cutoff_r,moved_atom)
        logical check_cutoffs
        real,dimension(:,:) ::cutoff_r
        integer moved_atom
        type(model) m
        integer, dimension(:), pointer :: atoms
        integer  nlist
        integer istat
        real radius, temp_x, temp_y, temp_z
        integer i,j
        integer num1, num2
        real dist_pair

        !find the maximum cut-off distance
        radius=maxval(maxval(cutoff_r, 1),1)

        ! Find all atoms within radius radius of moved_atom and put them 
        ! in the list atoms. Also sets nlist == size(atoms)+1.
        call hutch_list_3D(m, m%xx%ind(moved_atom),m%yy%ind(moved_atom),m%zz%ind(moved_atom), radius, atoms, istat, nlist)
        if (istat .eq. 1) then
            print *, 'memory allocation fails!'
            return
        else if(istat .eq. -1) then
            print *, 'no atom is found'
            check_cutoffs = .true. 
            return
        endif

        !begin to calculate pair distance
        !note, there are (nlist-1) atoms in total in 'atoms'

        !first, determine the type of moved_atom 
        do i=1, m%nelements
            if(m%znum%ind(moved_atom) .eq. m%atom_type(i)) then
                num1 = i
                exit
            endif
        enddo
        
        do i=1, (nlist-1)
            !do not count the pair distance to itself
            ! The below if statement should be irrelevant. moved_atom isn't put
            ! in atoms by hutch_list_3d as far as i know.
            check_cutoffs = .true.  !jwh 032409
            if (atoms(i) .ne. moved_atom) then
                do j=1, m%nelements
                    if(m%atom_type(j) .eq. m%znum%ind(atoms(i))) then
                        num2 = j
                        exit
                    endif
                enddo !j=1, m%nelements
              
                !endif  !032409 - jwh
               
                !calculate the atomic distance
                !compare with cutoff_r
                temp_x = m%xx%ind(moved_atom) - m%xx%ind(atoms(i))  !pbc added - jwh 04/14/2009
                temp_y = m%yy%ind(moved_atom) - m%yy%ind(atoms(i))
                temp_z = m%zz%ind(moved_atom) - m%zz%ind(atoms(i))
                temp_x = temp_x - m%lx*anint(temp_x/m%lx)
                temp_y = temp_y - m%ly*anint(temp_y/m%ly)
                temp_z = temp_z - m%lz*anint(temp_z/m%lz)
                  
                dist_pair = temp_x**2 + temp_y**2 + temp_z**2
                dist_pair = sqrt(dist_pair)
                   
                if (dist_pair  .lt. cutoff_r(num1, num2)) then
                    !write(*,*)"DEBUG", dist_pair  , cutoff_r(num1, num2) !debug - jwh 032409
                    check_cutoffs=.false.
                    exit
                endif
            endif !032409 - jwh
        enddo !i=1, (nlist-1)

        !if (associated(atoms)) deallocate(atoms) ! pmv 4/17/09
        if (nlist .gt. 1) deallocate(atoms) ! fy 4/17/09
    
    end function check_cutoffs
  


  !subroutine random_move---------------------------------------------------------------
  !Generates a random move of one atom.
  !Written for testing purpose here.
  !The mc move will be done some place else later, but it should still update the variables xx_cur, and xx_new (and etc).
  subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
    use mpi  
    type(model), intent(inout) :: m
    integer iseed, w
    real alpha, aa, bb, cc
    real, intent(out) :: xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new    !Cur and new positions of atoms
    real :: rand1, rand2, rand3, rand4


    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    rand3 = ran2(iseed)
    rand4 = ran2(iseed)

    w = int(m%natoms*rand1)+1

    !write(*,*)myid, rand1, rand2, rand3, rand4

    xx_cur = m%xx%ind(w)            !Cur positions of the atom before random move
    yy_cur = m%yy%ind(w)
    zz_cur = m%zz%ind(w)
    
    aa = alpha*(rand2 - 0.5)
    bb = alpha*(rand3 - 0.5)
    cc = alpha*(rand4 - 0.5)
    
    m%xx%ind(w) = m%xx%ind(w) + aa 
    m%yy%ind(w) = m%yy%ind(w) + bb 
    m%zz%ind(w) = m%zz%ind(w) + cc
    
    if(m%xx%ind(w)>m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)-m%lx       !pbc 
    if(m%yy%ind(w)>m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)-m%ly
    if(m%zz%ind(w)>m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)-m%lz
    if(m%xx%ind(w)<-m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)+m%lx
    if(m%yy%ind(w)<-m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)+m%ly
    if(m%zz%ind(w)<-m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)+m%lz
    
    xx_new=m%xx%ind(w)              !new positions of the atom after random move
    yy_new=m%yy%ind(w)
    zz_new=m%zz%ind(w)
  end subroutine random_move


END MODULE hrmc_functions
