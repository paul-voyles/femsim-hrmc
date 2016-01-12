! eam module: subroutines for calculating initial eam potential and updating energy after the mc move.
! Written by Jinwoo Hwang 08/10/2009
! Cleaned up by Jason on 07/02/2013

module eam_mod

    use model_mod 
    implicit none
    real, dimension(:,:), pointer :: f
    real, dimension(:,:), pointer :: rho
    real, dimension(:,:,:), pointer :: phi
    double precision, dimension(:), pointer :: e1, e2
    real :: drho, dr, eam_max_r
    logical, dimension(:), pointer :: not_counted
    integer :: nrho,  nr, nelements
                                  
contains

    ! Reads tabulated eam potential data
    ! Only called once in hrmc to initialize data.
    subroutine read_eam(m,eamfile)
        implicit none
        type(model), intent(in) :: m
        character (len=80) :: eamfile
        integer :: i, j, k, w, q, ii, jj
        integer, dimension(:), allocatable :: znum
        integer :: istat
        real :: rr, rho1
        real, dimension(:), allocatable :: mass, latt_const
        real, dimension(:,:,:), allocatable  :: f_temp
        real, dimension(:,:,:), allocatable  :: rho_temp
        real, dimension(:,:,:,:), allocatable :: phi_temp
        real, dimension(:,:), allocatable  :: f_temp2
        real, dimension(:,:), allocatable  :: rho_temp2
        real, dimension(:,:,:), allocatable :: phi_temp2
        real :: trash
        integer, dimension(:), allocatable :: reorder
        integer, dimension(:), allocatable :: atomic_numbers
        character (len=2), dimension(:), allocatable :: atom_syms
        character (len=2), dimension(103) :: all_atom_syms
        integer :: line

        all_atom_syms = (/'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

        open(unit=71,file=trim(eamfile),iostat=istat,form='formatted',status='unknown')
        read(71,*)  !comment line
        read(71,*)  !comment line
        read(71,*)  !comment line
        read(71,*) nelements
        allocate(atom_syms(nelements))
        allocate(atomic_numbers(nelements))
        allocate(reorder(nelements))
        rewind(71)
        read(71,*)  !comment line
        read(71,*)  !comment line
        read(71,*)  !comment line
        read(71,*) nelements, atom_syms
        read(71,*) nrho, drho, nr, dr, eam_max_r

        do i=1, nelements
            do j=1, 103
                if( trim(atom_syms(i)) .eq. trim(all_atom_syms(j)) ) then
                    atomic_numbers(i) = j
                endif
            enddo
        enddo
        do i=1, nelements
            do j=1, m%nelements
                if(atomic_numbers(i) .eq. m%atom_type(j)) then
                    reorder(i) = j
                endif
            enddo
        enddo

        if(myid.eq.0) write(*,*) "Atom types in EAM file: sym, Z, reorder#"
        do i=1, nelements
            if(myid.eq.0) write(*,*) "  ",atom_syms(i), atomic_numbers(i), reorder(i)
        enddo

        allocate(znum(nelements))
        allocate(mass(nelements))
        allocate(latt_const(nelements))
        allocate(f(nelements, nrho))
        allocate(rho(nelements, nr))
        allocate(phi(nelements, nelements, nr+1))

        ! WARNING! You probably cannot use the potential file directly
        ! that you downloaded. Rather, each number must be on a
        ! separate line. Use util/reformat_potential.py to convert.

        ! The elements in the eam file may not be in increasing
        ! atomic order (which is how znumr is set up) so we need
        ! to reorder if necessary. I set this up at the beginning.
        do i=1, nelements
            read(71,*) znum(i), mass(i), latt_const(i)
            do k=1, nrho
                read(71, *) f(reorder(i), k)
            enddo
            do k=1, nr
                read(71, *) rho(reorder(i), k)
            enddo
        enddo

        do i=1, nelements
            do j=1, nelements
                if(i.ge.j)then
                    do k=1, nr
                        read(71, *) phi(reorder(i), reorder(j), k)
                        phi(reorder(i), reorder(j), k) = phi(reorder(i), reorder(j), k)/(k*dr)
                        phi(reorder(j), reorder(i), k) = phi(reorder(i), reorder(j), k)
                    enddo
                endif
            enddo
        enddo
        close(71)

        allocate(e1(m%natoms))
        allocate(e2(m%natoms))
        allocate(not_counted(m%natoms))

        deallocate(znum)
        deallocate(mass, latt_const)
    end subroutine read_eam


    ! Calculates initial energy of model using eam potential
    ! Only called once in hrmc. Not used in Femsim of course.
    subroutine eam_initial(m, te1)
        implicit none
        type(model), intent(in) :: m
        integer :: i, j
        integer,  dimension(:), pointer :: atoms
        integer:: nlist, rbin, rhobin
        real :: xij, yij, zij, r, r2
        real :: phi1, phi2, rho1, rho2
        double precision, intent(out) :: te1
        integer :: istat
        do i=1, m%natoms
        phi2 = 0.0
        rho2 = 0.0
            call hutch_list_3D(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), eam_max_r, atoms, istat, nlist)
            do j=1, nlist-1
                if(atoms(j).ne.i)then
                    xij = m%xx%ind(i) - m%xx%ind(atoms(j))
                    yij = m%yy%ind(i) - m%yy%ind(atoms(j))          
                    zij = m%zz%ind(i) - m%zz%ind(atoms(j))
                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))
                    r2 = xij**2+yij**2+zij**2
                    if(r2.lt.eam_max_r*eam_max_r)then
                        r = sqrt(r2)
                        rbin = int(r/dr)+1
                        if(rbin == nr+1) rbin = rbin - 1 ! Jason due to out of bounds error. I believe this is the correct fix, and that the problem arrises due to a rounding error. 20130731
                        phi1 = phi(m%znum_r%ind(i), m%znum_r%ind(atoms(j)) , rbin)
                        phi2 = phi2 + phi1
                        rho1 = rho(m%znum_r%ind(i), rbin)
                        rho2 = rho2 + rho1
                    endif
                endif
            enddo
            rhobin = int(rho2/drho)+1
            if(rhobin.le.0)then
                e1(i) = 0.5*phi2
            else
                e1(i) = f(m%znum_r%ind(i), rhobin) + 0.5*phi2
            endif
            if(associated(atoms)) deallocate(atoms)
        enddo

        te1=0.0
        do i=1, m%natoms
            te1 = te1 + e1(i)
        enddo
        !write(*,*)"initial energy=", te1
    end subroutine eam_initial


    ! Calculates initial energy of model using eam potential without hutch
    ! Currently never called.
    subroutine eam_initial_no_hutch(m, te1)
        implicit none
        type(model), intent(in) :: m
        integer :: i, j
        integer:: nlist, rbin, rhobin
        real :: xij, yij, zij, r, r2
        real :: phi1, phi2, rho1, rho2
        real, intent(out) :: te1
        integer :: istat

        allocate(e1(m%natoms))
        allocate(e2(m%natoms))
        do i=1, m%natoms
            phi2 = 0.0
            rho2 = 0.0
            do j=1, m%natoms
                if(j.ne.i)then
                    xij = m%xx%ind(i) - m%xx%ind(j)
                    yij = m%yy%ind(i) - m%yy%ind(j)         
                    zij = m%zz%ind(i) - m%zz%ind(j)
                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))
                    r2 = xij**2+yij**2+zij**2
                    if(r2.le.eam_max_r*eam_max_r)then
                        r = sqrt(r2)
                        rbin = int(r/dr)+1
                        if(rbin == nr+1) rbin = rbin - 1 ! Jason due to out of bounds error. I believe this is the correct fix, and that the problem arrises due to a rounding error. 20130731
                        phi1 = phi(m%znum_r%ind(i), m%znum_r%ind(j) , rbin)
                        phi2 = phi2 + phi1
                        rho1 = rho(m%znum_r%ind(i), rbin)
                        rho2 = rho2 + rho1
                    endif
                endif
            enddo
            rhobin = int(rho2/drho)+1
            if(rhobin.le.0)then
                e1(i) = 0.5*phi2
            else
                e1(i) = f(m%znum_r%ind(i), rhobin) + 0.5*phi2
            endif
        enddo

        te1=0.0
        do i=1, m%natoms
            te1 = te1 + e1(i)
        enddo
        !write(*,*)"initial energy=", te1
    end subroutine eam_initial_no_hutch


    ! Calculates initial energy of model using eam potential
    ! Called in hrmc loop.
    subroutine eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
        implicit none
        type(model), intent(in) :: m
        integer :: i, j, k
        integer, intent(in) :: w
        integer,  dimension(:), pointer :: atoms1, atoms2, atoms3, atoms4
        integer:: nlist1, nlist2, nlist3, nlist4, rbin, rhobin
        real :: xij, yij, zij, r, r2
        real :: phi1, phi2, rho1, rho2
        double precision, intent(out) :: te2
        real, intent(in) :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new 
        integer :: istat
        not_counted = .true.

        call hutch_list_3D(m, xx_cur, yy_cur, zz_cur, eam_max_r, atoms1, istat, nlist1)

        do i=1, nlist1-1
            phi2 = 0.0
            rho2 = 0.0
            call hutch_list_3D(m, m%xx%ind(atoms1(i)), m%yy%ind(atoms1(i)), m%zz%ind(atoms1(i)), eam_max_r, atoms2, istat, nlist2)
            do j=1, nlist2-1
                if(atoms2(j).ne.atoms1(i))then
                    xij = m%xx%ind(atoms1(i)) - m%xx%ind(atoms2(j))
                    yij = m%yy%ind(atoms1(i)) - m%yy%ind(atoms2(j))         
                    zij = m%zz%ind(atoms1(i)) - m%zz%ind(atoms2(j))
                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))
                    r2 = xij**2+yij**2+zij**2
                    if(r2.le.eam_max_r*eam_max_r)then
                        r = sqrt(r2)
                        rbin = int(r/dr)+1
                        if(rbin == nr+1) rbin = rbin - 1 ! Jason due to out of bounds error. I believe this is the correct fix, and that the problem arrises due to a rounding error. 20130731
                        phi1 = phi(m%znum_r%ind(atoms1(i)), m%znum_r%ind(atoms2(j)) , rbin)
                        phi2 = phi2 + phi1
                        rho1 = rho(m%znum_r%ind(atoms1(i)), rbin)
                        rho2 = rho2 + rho1
                    endif
                endif
            enddo
            rhobin = int(rho2/drho)+1
            if(rhobin.le.0)then
                e2(atoms1(i)) = 0.5*phi2
            else
                e2(atoms1(i)) = f(m%znum_r%ind(atoms1(i)), rhobin) + 0.5*phi2
            endif
            not_counted(atoms1(i)) = .false.
            if(associated(atoms2)) deallocate(atoms2)
        enddo

        call hutch_list_3D(m, xx_new, yy_new, zz_new, eam_max_r, atoms3, istat, nlist3)

        do i=1, nlist3-1
            phi2 = 0.0
            rho2 = 0.0
            if(not_counted(atoms3(i)))then
                call hutch_list_3D(m, m%xx%ind(atoms3(i)), m%yy%ind(atoms3(i)), m%zz%ind(atoms3(i)), eam_max_r, atoms4, istat, nlist4)
                do j=1, nlist4-1
                    if(atoms4(j).ne.atoms3(i))then
                        xij = m%xx%ind(atoms3(i)) - m%xx%ind(atoms4(j))
                        yij = m%yy%ind(atoms3(i)) - m%yy%ind(atoms4(j))     
                        zij = m%zz%ind(atoms3(i)) - m%zz%ind(atoms4(j))
                        xij = xij-m%lx*anint(xij/(m%lx))                
                        yij = yij-m%ly*anint(yij/(m%ly))            
                        zij = zij-m%lz*anint(zij/(m%lz))
                        r2 = xij**2+yij**2+zij**2
                        if(r2.le.eam_max_r*eam_max_r)then
                            r = sqrt(r2)
                            rbin = int(r/dr)+1
                            if(rbin == nr+1) rbin = rbin - 1 ! Jason due to out of bounds error. I believe this is the correct fix, and that the problem arrises due to a rounding error. 20130731
                            phi1 = phi(m%znum_r%ind(atoms3(i)), m%znum_r%ind(atoms4(j)) , rbin)
                            phi2 = phi2 + phi1
                            rho1 = rho(m%znum_r%ind(atoms3(i)), rbin)
                            rho2 = rho2 + rho1
                        endif
                    endif
                enddo
                rhobin = int(rho2/drho)+1
                if(rhobin.le.0)then
                    e2(atoms3(i)) = 0.5*phi2
                else
                    e2(atoms3(i)) = f(m%znum_r%ind(atoms3(i)), rhobin) + 0.5*phi2
                endif
                if(associated(atoms4)) deallocate(atoms4)
            endif
        enddo

        te2=0.0
        do i=1, m%natoms
            te2 = te2 + e2(i)
        enddo
        if(associated(atoms1)) deallocate(atoms1)
        if(associated(atoms3)) deallocate(atoms3)
    end subroutine eam_mc

end module eam_mod
