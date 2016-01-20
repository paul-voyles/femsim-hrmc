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

    function chi_square(alpha,vk_exp, vk_exp_err, vk, scale_fac, nk)
        real :: alpha
        double precision, pointer, dimension(:) :: vk_exp, vk_exp_err
        double precision, pointer, dimension(:) :: vk
        real, intent(in) :: scale_fac
        integer, intent(in) ::nk
        integer j
        integer nf   ! normalization factor
        double precision :: chi_square
        chi_square=0.0
        
        nf = 0
        do j=1,nk
            nf = nf + 1
            chi_square = chi_square+alpha*((vk_exp(j)*scale_fac-vk(j))/(scale_fac*vk_exp_err(j)))**2
        enddo
        chi_square = chi_square/nf
    end function chi_square


    ! Returns .true. if no two atoms are within cutoff distance.
    ! Otherwise returns .false.
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

        ! Find the maximum cut-off distance
        radius = maxval(maxval(cutoff_r, 1),1)

        ! Find all atoms within radius radius of moved_atom and put them 
        ! in the list atoms. Also sets nlist == size(atoms)+1.
        call hutch_list_3D(m, m%xx%ind(moved_atom), m%yy%ind(moved_atom), m%zz%ind(moved_atom), radius, atoms, istat, nlist)
        if (istat .eq. 1) then
            write(0, *) 'memory allocation fails!'
            return
        else if(istat .eq. -1) then
            write(0, *) 'no atom is found'
            check_cutoffs = .true. 
            return
        endif

        ! First, determine the type of moved_atom so we can use it in cutoff_r
        do i=1, m%nelements
            if(m%znum%ind(moved_atom) .eq. m%atom_type(i)) then
                num1 = i
                exit
            endif
        enddo
        
        do i=1, (nlist-1)
            ! Do not count the pair distance to itself
            ! The below if statement should be irrelevant. moved_atom isn't put
            ! in atoms by hutch_list_3d as far as i know. However, it's a nice
            ! safety feature and pretty irrelelvant.
            check_cutoffs = .true.
            if (atoms(i) .ne. moved_atom) then
                do j=1, m%nelements
                    if(m%atom_type(j) .eq. m%znum%ind(atoms(i))) then
                        num2 = j
                        exit
                    endif
                enddo !j=1, m%nelements
              
                !calculate the atomic distance and compare with cutoff_r
                temp_x = m%xx%ind(moved_atom) - m%xx%ind(atoms(i))  ! periodic boundary conditions
                temp_y = m%yy%ind(moved_atom) - m%yy%ind(atoms(i))
                temp_z = m%zz%ind(moved_atom) - m%zz%ind(atoms(i))
                temp_x = temp_x - m%lx*anint(temp_x/m%lx)
                temp_y = temp_y - m%ly*anint(temp_y/m%ly)
                temp_z = temp_z - m%lz*anint(temp_z/m%lz)
                  
                dist_pair = temp_x**2 + temp_y**2 + temp_z**2
                dist_pair = sqrt(dist_pair)
                   
                if (dist_pair  .lt. cutoff_r(num1, num2)) then
                    check_cutoffs=.false.
                    exit
                endif
            endif
        enddo

        if (nlist .gt. 1) deallocate(atoms)
    end function check_cutoffs
  


  ! Generates a random move of one atom.
  ! The actual positions are not changed, just returned.
  subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
    use mpi  
    type(model), intent(inout) :: m
    integer iseed, w
    real alpha, aa, bb, cc
    real, intent(out) :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: rand1, rand2, rand3, rand4


    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    rand3 = ran2(iseed)
    rand4 = ran2(iseed)

    w = int(m%natoms*rand1)+1

    ! Current positions of the atom before random move
    xx_cur = m%xx%ind(w)
    yy_cur = m%yy%ind(w)
    zz_cur = m%zz%ind(w)
    
    aa = alpha*(rand2 - 0.5)
    bb = alpha*(rand3 - 0.5)
    cc = alpha*(rand4 - 0.5)
    
    m%xx%ind(w) = m%xx%ind(w) + aa 
    m%yy%ind(w) = m%yy%ind(w) + bb 
    m%zz%ind(w) = m%zz%ind(w) + cc
    
    ! Periodic boundary conditions
    if(m%xx%ind(w) >  m%lx*0.5) m%xx%ind(w) = m%xx%ind(w) - m%lx
    if(m%yy%ind(w) >  m%ly*0.5) m%yy%ind(w) = m%yy%ind(w) - m%ly
    if(m%zz%ind(w) >  m%lz*0.5) m%zz%ind(w) = m%zz%ind(w) - m%lz
    if(m%xx%ind(w) < -m%lx*0.5) m%xx%ind(w) = m%xx%ind(w) + m%lx
    if(m%yy%ind(w) < -m%ly*0.5) m%yy%ind(w) = m%yy%ind(w) + m%ly
    if(m%zz%ind(w) < -m%lz*0.5) m%zz%ind(w) = m%zz%ind(w) + m%lz
    
    ! New positions of the atom after random move
    xx_new = m%xx%ind(w)
    yy_new = m%yy%ind(w)
    zz_new = m%zz%ind(w)
  end subroutine random_move


END MODULE hrmc_functions
