!Gr module: subroutines for calculating initial gr and updating gr after the mc move.

! Subroutines:

!subroutine scatt_power
!Calculates electron, neutron, xray scattering power depending on which experimental data were used.
!Uses 'module scattering factors'.

!subroutine gr_initialize
!Initialize partial and total grs of any kind if used. 
!Calculates radius from max range of experimental data.

!subroutine gr_no_hutch
!Calculates gr without hutch. Can be used for initial gr calculation and etc.

!subroutine gr_hutch_mc
!Updates gr with hutch with a moved atom in one mc step.
!Calls subroutine hutch_list_3D,
!Subtracts the pairs from the ole atom position, and after the mc move, adds the pairs from the new atom position.

!subroutine accept_gr
!Updates gr with new ones.
!Should be called when hrmc accepts the new position.

!subroutine reject_gr
!Returns gr to the old one.
!Should be called when hrmc rejects the new position.

!subroutine gr_hutch
!Does the same thing as the subroutine gr_no_hutch, but this subroutine uses hutch algorithm.
!Written for testing purpose here.
!Calculates gr with hutch.


! Changelog
!
! Gr module:
! Written by Jinwoo Hwang 01/09/2009

! Modified by Jinwoo Hwang 01/13/2009

! Subroutines scatt_power, gr_initialize, gr_no_hutch, gr_hutch, gr_hutch_mc, accept_gr, 
! and reject_gr were modified to contain all 3 types of gr - electron, neutron, and xray. 
! They now use logical parameter"used_data_set", which is defined in the main program, 
! to determine which types of gr to calculate. The mbin (# of bins in gr) is devided to mbin_e, mbin_n, and mbin_x
! to count on the difference in number of bins among the experimental data. 
! Same thing done for the bin size, del_r. 
! The xray scattering power is calculated using xray scattering factor. - Jinwoo Hwang 01/18/2009
! 
! subroutines used to calculate electron scatter factors 
! and X-ray scatter factors
! take zero or fewer arguments than previous version
! updated by Feng Yi on 01/25/2009
! subroutine reduced_znum moved to model_module - Jinwoo Hwang 01/27/2009

!revised by Feng Yi on 02/27/2009, check comment with 02/27/2009
!
!
!subroutine gr_hutch_mc revided by Feng Yi on 03/03/2009 
!
!subroutine gr_hutch_mc revised by Feng Yi on 03/04/2009
!
!dev_del_r is added to the gr_no_hutch JWH 03/30/2009
!**********************************************************
module gr_mod

    use model_mod 

    implicit none


    integer mbin_e,mbin_n,mbin_x        !mbin=number of bins in gr. w=random number                         

    real, pointer, dimension(:,:) :: be,bn,bx   !Scattering power for each pair type.(electron,neutron,xray)
    INTEGER, pointer, dimension(:,:,:) :: gr_e_part_cur,gr_e_part_new  !Simulated partial pdf
    INTEGER, pointer, dimension(:,:,:) :: gr_n_part_cur,gr_n_part_new  !Simulated partial pdf
    INTEGER, pointer, dimension(:,:,:) :: gr_x_part_cur,gr_x_part_new  !Simulated partial pdf
    real, pointer, dimension(:) :: gr_e_sim_cur,gr_e_sim_new        !Simulated total e gr
    real, pointer, dimension(:) :: gr_n_sim_cur,gr_n_sim_new        !Simulated total n gr
    real, pointer, dimension(:) :: gr_x_sim_cur,gr_x_sim_new        !Simulated total x gr
    real del_r_e, del_r_n, del_r_x     !Bin size in gr. 
    real radius
                                      
contains


    !subroutine scatt_power
    !Calculates electron, neutron, xray scattering power depending on which experimental data were used.
    !Uses 'module scattering factors'.
    subroutine scatt_power(m,used_data_sets,istat)
        use scattering_factors
        type(model), intent(in) :: m                                
        integer, intent(out) :: istat
        integer i, j, k                        
        !real, dimension(103,3):: a_e, b_e, c_e, d_e
        !real, dimension(98,4):: a_x, b_x, c_x, d_x
        real q
        real scatfact_e(103,100),scatfact_x(98,100)
        real scat_binsize
        real, allocatable, dimension(:,:,:) :: fifj
        real, allocatable, dimension(:,:,:) :: b_q     !Scattering power as a function of q.    
        real, dimension(100) :: xifi
        logical, dimension(4), intent(in) :: used_data_sets

        allocate(fifj(m%nelements,m%nelements,100), stat=istat)
        allocate(b_q(m%nelements,m%nelements,100), stat=istat)

        if(used_data_sets(1)) then   
            allocate(be(m%nelements,m%nelements), stat=istat)

            scat_binsize = 0.04
            
            !call read_f_e(a_e, b_e, c_e, d_e)
            call read_f_e ! because this subroutine does not take argument after update
                          ! by fy on 01/25/2009

            ! This is a 2D array. It is defined a few lines above.
            scatfact_e = 0.0
     
            do j=1,m%nelements
                do k=1, 100
                        q=k*scat_binsize
                        !scatfact_e(j,k)=f_e(m%atom_type(j),q,a_e, b_e, c_e, d_e)
                        scatfact_e(j,k)=f_e(m%atom_type(j),q) !updated by fy on 01/25/2009
                enddo
            enddo
        
            xifi = 0.0

            do k=1, 100
                do i=1, m%nelements
                    xifi(k) = xifi(k) + m%composition(i)*scatfact_e(i,k)
                enddo
            enddo
            
            fifj = 0.0
            b_q = 0.0
            be = 0.0
            
            do k=1, 100
                do i=1, m%nelements
                    do j=1, m%nelements
                        fifj(i,j,k)=scatfact_e(i,k)*scatfact_e(j,k)
                        b_q(i,j,k)=fifj(i,j,k)/(xifi(k)**2)
                        be(i,j)=be(i,j) + b_q(i,j,k)/100.0
                    enddo
                enddo
            enddo
        
        endif

        if(used_data_sets(2)) then

            bn = 0.0
            allocate(bn(m%nelements,m%nelements), stat=istat)
            do i=1, m%nelements
                do j=1, m%nelements
                    bn(i,j)=1.0         !temporary
                enddo
            enddo
        endif

        if(used_data_sets(3)) then

            allocate(bx(m%nelements,m%nelements), stat=istat)

            scat_binsize = 0.04
            
            !call read_f_x(a_x, b_x, c_x)
            call read_f_x !updated by fy on 01/25/2009

            scatfact_x = 0.0
        
            do j=1,m%nelements     
                do i=1, 4
                    do k=1, 100
                        q=k*scat_binsize
                        !scatfact_x(j,k)=f_x(m%atom_type(j),q,a_x, b_x, c_x)
                        scatfact_x(j,k)=f_x(m%atom_type(j),q) !updated by fy on 01/25/2009
                    enddo
                enddo
            enddo
        
            xifi = 0.0
            do k=1, 100
                do i=1, m%nelements
                    xifi(k) = xifi(k) + m%composition(i)*scatfact_x(i,k)
                enddo
            enddo
            
            fifj = 0.0
            b_q = 0.0
            bx = 0.0
            
            do k=1, 100
                do i=1, m%nelements
                    do j=1, m%nelements
                        fifj(i,j,k)=scatfact_x(i,k)*scatfact_x(j,k)
                        b_q(i,j,k)=fifj(i,j,k)/(xifi(k)**2)
                        bx(i,j)=bx(i,j) + b_q(i,j,k)/100.0
                    enddo
                enddo
            enddo
        endif
        deallocate(fifj,b_q)

    end subroutine scatt_power



    !subroutine gr_initialize
    !Initialize partial and total grs of any kind if used. 
    !Calculates radius from max range of experimental data.
    subroutine gr_initialize(m,r_e,gr_e,r_n,gr_n,r_x,gr_x,used_data_sets,istat)
        implicit none   
        type(model), intent(in) :: m
        integer, intent(out) :: istat
        integer i,j,k
        real, dimension(:), intent(in) :: r_e, gr_e, r_n, gr_n, r_x, gr_x
        real start_r_e, start_r_n, start_r_x
        logical, dimension(4), intent(in) :: used_data_sets
        integer mbin_e1, mbin_n1, mbin_x1 !added by Feng Yi on 02/27/2009

        if(used_data_sets(1)) then
            mbin_e = size(r_e) 
            mbin_e1 = mbin_e + 1
            start_r_e = r_e(1)
            del_r_e = r_e(2) - r_e(1)
        
            ! Allocate two arrays.
            allocate(gr_e_sim_cur(mbin_e1), stat=istat) !02/27/2009
            allocate(gr_e_sim_new(mbin_e1), stat=istat) !02/27/2009

            gr_e_sim_cur = 0.0
            gr_e_sim_new = 0.0

            allocate(gr_e_part_cur(m%nelements, m%nelements,mbin_e1), stat=istat) !modified by Feng on 02/27/2009
            allocate(gr_e_part_new(m%nelements, m%nelements,mbin_e1), stat=istat) !02/27/2009
        
            do i=1,mbin_e1 !02/27/2009
                do j=1,m%nelements
                    do k=1,m%nelements
                        gr_e_part_cur(j,k,i) = 0
                        gr_e_part_new(j,k,i) = 0
                    enddo
                enddo
            enddo
        endif

        if(used_data_sets(2)) then
            mbin_n = size(r_n) 
            mbin_n1 = mbin_n + 1 !added by Feng Yi on 02/27/2009
            start_r_n = r_n(1)
            del_r_n = r_n(2) - r_n(1)
            
        
            allocate(gr_n_sim_cur(mbin_n1), stat=istat) !02/27/2009
            allocate(gr_n_sim_new(mbin_n1), stat=istat) !02/27/2009

            gr_n_sim_cur = 0.0
            gr_n_sim_new = 0.0

            allocate(gr_n_part_cur(m%nelements, m%nelements,mbin_n1), stat=istat) !02/27/2009
            allocate(gr_n_part_new(m%nelements, m%nelements,mbin_n1), stat=istat) !02/27/2009
        
            do i=1,mbin_n1 !02/27/2009
                do j=1,m%nelements
                    do k=1,m%nelements
                        gr_n_part_cur(j,k,i) = 0.0
                        gr_n_part_new(j,k,i) = 0.0
                    enddo
                enddo
            enddo
        endif

        if(used_data_sets(3)) then
            mbin_x = size(r_x) + 1 
            !write(*,*) mbin_x
            mbin_x1 = mbin_x + 1 !added by Feng Yi on 02/27/2009
            start_r_x = r_x(1)
            del_r_x = r_x(2) - r_x(1)
            
        
            allocate(gr_x_sim_cur(mbin_x1), stat=istat) !02/27/2009
            allocate(gr_x_sim_new(mbin_x1), stat=istat) !02/27/2009

            gr_x_sim_cur = 0.0
            gr_x_sim_new = 0.0

            allocate(gr_x_part_cur(m%nelements, m%nelements,mbin_x1), stat=istat) !02/27/2009
            allocate(gr_x_part_new(m%nelements, m%nelements,mbin_x1), stat=istat) !02/27/2009
        
            do i=1,mbin_x1 !02/27/2009
                do j=1,m%nelements
                    do k=1,m%nelements
                        gr_x_part_cur(j,k,i) = 0.0
                        gr_x_part_new(j,k,i) = 0.0
                    enddo
                enddo
            enddo
        endif

        radius = max(mbin_e*del_r_e,mbin_n*del_r_n,mbin_x*del_r_x)
    end subroutine gr_initialize



    !subroutine gr_no_hutch
    !Calculates gr without hutch. Can be used for initial gr calculation and etc.
    subroutine gr_no_hutch(m,used_data_sets)
        type(model), intent(in) :: m
        integer i, j, k, ig
        real rmax_e,rmax_n,rmax_x
        real vb, nid, rho
        real xi, yi, zi, xj, yj, zj, xij, yij, zij, R, R2
        real dev_del_r_e, dev_del_r_n, dev_del_r_x
        logical, dimension(4), intent(in) :: used_data_sets

        rho = m%natoms/(m%lx*m%ly*m%lz)  

        dev_del_r_e = del_r_e / 100.0
        dev_del_r_n = del_r_n / 100.0
        dev_del_r_x = del_r_x / 100.0                       

        do i=1, m%natoms
            !WRITE(*,*) 'i is: ', i !debug
            xi = m%xx%ind(i)
            yi = m%yy%ind(i)
            zi = m%zz%ind(i)

            do j=1, m%natoms
                if(i.ne.j)then
                    !WRITE(*,*) 'j is: ', j !debug
                    xj = m%xx%ind(j)
                    yj = m%yy%ind(j)
                    zj = m%zz%ind(j)
                
                    xij = xi - xj
                    yij = yi - yj
                    zij = zi - zj

                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))            
        
                    R2 = xij**2+yij**2+zij**2
                    R = sqrt(R2)

                    if(used_data_sets(1)) then
                        rmax_e = del_r_e*mbin_e
                        if(R.lt.rmax_e)then
                            ig = int(R/del_r_e) +1
                            IF(ig .GT. (mbin_e )) THEN
                                write(*,*) 'Initialization error!', ig, mbin_e 
                            endif
                            !WRITE(*,*) 'Wrong here!' !debug
                            !WRITE(*,*) 'ig is: ', ig, 'mbin_e is; ', mbin_e
                            gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)+1
                            gr_e_part_new(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)
                        endif
                    endif

                    if(used_data_sets(2)) then
                        rmax_n = del_r_n*mbin_n
                        if(R.lt.rmax_n)then
                            ig = int(R/del_r_n) +1
                            gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)+1
                            gr_n_part_new(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)
                        endif
                    endif

                    if(used_data_sets(3)) then
                        rmax_x = del_r_x*mbin_x
                        if(R.lt.rmax_x)then
                            ig = int(R/del_r_x) +1
                            gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)+1
                            gr_x_part_new(m%znum_r%ind(i),m%znum_r%ind(j),ig) = gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(j),ig)
                        endif
                    endif
                endif
            enddo
        enddo

        if(used_data_sets(1)) then
            do i=1,mbin_e
                R = del_r_e*(i-1.0) + dev_del_r_e
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_e_sim_cur(i) = gr_e_sim_cur(i) + be(j,k)*gr_e_part_cur(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=903,file="test_gr_initial.txt",form='formatted',status='unknown')
            do i=1, mbin_e
                R = del_r_e*(i)-del_r_e + dev_del_r_e
                vb =((i+1)**3-i**3)*del_r_e**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_e_sim_cur(i) = (gr_e_sim_cur(i)/(m%natoms)-nid)/(R*del_r_e)                
                !write(903,*)R, gr_e_sim_cur(i)
            enddo
            !close(903)
        endif

        if(used_data_sets(2)) then
            do i=1,mbin_n
                R = del_r_n*(i-1.0) + dev_del_r_n
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_n_sim_cur(i) = gr_n_sim_cur(i) + bn(j,k)*gr_n_part_cur(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=904,file="test_tot_n_gr.txt",form='formatted',status='unknown')
            do i=1, mbin_n
                R = del_r_n*(i)-del_r_n + dev_del_r_n
                vb =((i+1)**3-i**3)*del_r_n**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_n_sim_cur(i) = (gr_n_sim_cur(i)/(m%natoms)-nid)/(R*del_r_n)               
                !write(904,*)R, gr_n_sim_cur(i)
            enddo
            !close(904)
        endif

        if(used_data_sets(3)) then
            do i=1,mbin_x
                R = del_r_x*(i-1.0) + dev_del_r_x
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_x_sim_cur(i) = gr_x_sim_cur(i) + bx(j,k)*gr_x_part_cur(j,k,i) 
                    enddo
                enddo
            enddo
            !write(*,*)"bx=", bx(1,1)
            !write(*,*)"mbin_x=", mbin_x
            !write(*,*)"R", R
            !open(unit=905,file="test_tot_x_gr.txt",form='formatted',status='unknown')
            do i=1, mbin_x
                R = del_r_x*(i)-del_r_x + dev_del_r_x
                vb =((i+1)**3-i**3)*del_r_x**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_x_sim_cur(i) = (gr_x_sim_cur(i)/(m%natoms)-nid)/(R*del_r_x)                
            !   write(905,*)R, gr_x_sim_cur(i)
            enddo
            !close(905)
        endif
    end subroutine gr_no_hutch


    !subroutine gr_hutch_mc
    !Updates gr with hutch with a moved atom in one mc step.
    !Calls subroutine hutch_list_3D,
    !Subtracts the pairs from the ole atom position, and after the mc move, adds the pairs from the new atom position. 
    subroutine gr_hutch_mc(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new,used_data_sets,istat)   
        type(model), intent(inout) :: m
        real vb, nid, rho
        real xi, yi, zi, xj, yj, zj, xij, yij, zij, R, R2
        ! Current and new  positions of atoms, by Feng Yi on 03/03/2009
        real, INTENT(IN) :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
        integer, intent(inout) :: istat
        integer,  dimension(:), pointer :: atoms
        integer:: nlist 
        integer i, j, k, ig ,w
        logical, dimension(4), intent(in) :: used_data_sets
        REAL dev_del_r_e, dev_del_r_n, dev_del_r_x !added by Feng Yi on 02/27/2009
        !solve divided by 0 problem

        rho = m%natoms/(m%lx*m%ly*m%lz)

        dev_del_r_e = del_r_e / 100.0
        dev_del_r_n = del_r_n / 100.0
        dev_del_r_x = del_r_x / 100.0

        m%xx%ind(w) = xx_cur
        m%yy%ind(w) = yy_cur
        m%zz%ind(w) = zz_cur !revisded by Feng Yi on 03/03/2009
        !uncommented by Feng Yi on 03/06/2009 for compatibility with random_move in hrmc_functions.f90

        call hutch_list_3D(m, m%xx%ind(w), m%yy%ind(w), m%zz%ind(w), radius, atoms, istat, nlist)

        do k=1,nlist-1
            if(atoms(k).ne.w)then

                xj = m%xx%ind(atoms(k))
                yj = m%yy%ind(atoms(k))
                zj = m%zz%ind(atoms(k))
            
                xij = m%xx%ind(w) - xj
                yij = m%yy%ind(w) - yj
                zij = m%zz%ind(w) - zj
            
                xij = xij-m%lx*anint(xij/(m%lx))                
                yij = yij-m%ly*anint(yij/(m%ly))            
                zij = zij-m%lz*anint(zij/(m%lz))            
                
                R2 = xij**2+yij**2+zij**2
                R=sqrt(R2)

                if(R.le.radius)then
                    if(used_data_sets(1)) then
                        ig = int(R/del_r_e) +1 !should + 1 NOT 1.0 Changed by Feng Yi on 02/27/2009
                        gr_e_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_e_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)-1   !revised by Feng Yi on 03/04/2009
                        gr_e_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_e_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)-1   !revised by Feng Yi on 03/05/2009
                    endif

                    if(used_data_sets(2)) then
                        ig = int(R/del_r_n) +1
                        gr_n_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_n_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)-1  !changed by Feng on 03/04/3009
                        gr_n_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_n_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)-1   !revised by Feng Yi on 03/05/2009
                    endif

                    if(used_data_sets(3)) then
                        ig = int(R/del_r_x) +1
                        gr_x_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_x_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)-1  !03/04/2009 by Feng Yi
                        gr_x_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_x_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)-1   !revised by Feng Yi on 03/05/2009
                    endif
                endif !R < radius
            endif !w != atoms(k)
        enddo

        !update position here
        m%xx%ind(w) = xx_new
        m%yy%ind(w) = yy_new
        m%zz%ind(w) = zz_new
        if(associated(atoms)) deallocate(atoms) !added by Feng Yi on 03/04/2009

        call hutch_list_3D(m, m%xx%ind(w), m%yy%ind(w), m%zz%ind(w), radius, atoms, istat, nlist)

        do k=1,nlist-1
            if(atoms(k).ne.w)then

                xj = m%xx%ind(atoms(k))
                yj = m%yy%ind(atoms(k))
                zj = m%zz%ind(atoms(k))
            
                xij = m%xx%ind(w) - xj
                yij = m%yy%ind(w) - yj
                zij = m%zz%ind(w) - zj
            
                !original code
                xij = xij-m%lx*anint(xij/(m%lx))                
                yij = yij-m%ly*anint(yij/(m%ly))            
                zij = zij-m%lz*anint(zij/(m%lz))            
                
                R2 = xij**2+yij**2+zij**2
                
                !code added by Feng Yi on 02/27/2009
                !xij = MIN((xij)**2, (xij-m%lx)**2, (xij+m%lx)**2)
                !yij = MIN((yij)**2, (yij-m%ly)**2, (yij+m%ly)**2)
                !zij = MIN((zij)**2, (zij-m%lz)**2, (zij+m%lz)**2)
                !R2 = xij + yij + zij
                !*******************************************

                R=sqrt(R2)

                if(R.le.radius)then
                    if(used_data_sets(1)) then
                        ig = int(R/del_r_e) +1
                        !IF(ig .GT. mbin_e) THEN
                        !  WRITE(*, *) 'outside mbin_e! ', ig, mbin_e
                        !ENDIF
                        gr_e_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_e_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)+1  !03/04/2009 by Feng Yi
                        gr_e_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_e_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)+1   !revised by Feng Yi on 03/05/2009
                    endif

                    if(used_data_sets(2)) then
                        ig = int(R/del_r_n) +1
                        gr_n_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_n_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)+1  !by Feng Yi on 03/04/2009
                        gr_n_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_n_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)+1   !revised by Feng Yi on 03/05/2009
                    endif

                    if(used_data_sets(3)) then
                        ig = int(R/del_r_x) +1
                        gr_x_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig) = gr_x_part_new(m%znum_r%ind(w),m%znum_r%ind(atoms(k)),ig)+1  !by feng yi on 03/04/2009
                        gr_x_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig) = gr_x_part_new(m%znum_r%ind(atoms(k)),m%znum_r%ind(w),ig)+1   !revised by Feng Yi on 03/05/2009
                    endif
                endif

            endif
        enddo

        !open(unit=907,file="test_gr_part_new.txt",form='formatted',status='unknown')
        !do i=1,mbin_e
        !   R = del_r_e*(i-1.0)
        !   write(906,*) r,gr_e_part_cur(1,1,i)
        !   write(907,*) r,gr_e_part_new(1,1,i)
        !enddo
        !close(906)
        !close(907)
        
        if(used_data_sets(1)) then
            !Initialize gr_e_sim_new by Feng Yi on 03/05/2009
            gr_e_sim_new = 0.0
            do i=1,mbin_e
                R = del_r_e*(i-1.0) + dev_del_r_e
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_e_sim_new(i) = gr_e_sim_new(i) + be(j,k)*gr_e_part_new(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=908,file="test_tot_e_gr_mc.txt",form='formatted',status='unknown')
            do i=1, mbin_e
                R = del_r_e*(i)-del_r_e + dev_del_r_e
                vb =((i+1)**3-i**3)*del_r_e**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_e_sim_new(i) = (gr_e_sim_new(i)/(m%natoms)-nid)/(R*del_r_e)                
                !write(908,*)R, gr_e_sim_new(i)
            enddo
            !close(908)
        endif

        if(used_data_sets(2)) then
            !Initialize gr_e_sim_new by Feng Yi on 03/05/2009
            gr_n_sim_new = 0.0
            do i=1,mbin_n
                R = del_r_n*(i-1.0) + dev_del_r_n
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_n_sim_new(i) = gr_n_sim_new(i) + bn(j,k)*gr_n_part_new(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=909,file="test_tot_n_gr_mc.txt",form='formatted',status='unknown')
            do i=1, mbin_n
                R = del_r_n*(i)-del_r_n + dev_del_r_n
                vb =((i+1)**3-i**3)*del_r_n**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_n_sim_new(i) = (gr_n_sim_new(i)/(m%natoms)-nid)/(R*del_r_n)               
                !write(909,*)R, gr_n_sim_new(i)
            enddo
            !close(909)
        endif


        if(used_data_sets(3)) then  
            !Initialize gr_e_sim_new by Feng Yi on 03/05/2009
            gr_x_sim_new = 0.0
            do i=1,mbin_x
                R = del_r_x*(i-1.0) + dev_del_r_x
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_x_sim_new(i) = gr_x_sim_new(i) + bx(j,k)*gr_x_part_new(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=911,file="test_tot_x_gr_mc.txt",form='formatted',status='unknown')
            do i=1, mbin_x
                R = del_r_x*(i)-del_r_x + dev_del_r_x
                vb =((i+1)**3-i**3)*del_r_x**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_x_sim_new(i) = (gr_x_sim_new(i)/(m%natoms)-nid)/(R*del_r_x)                
                !write(911,*)R, gr_x_sim_new(i)
            enddo
            !close(911)
        endif
        
        if(associated(atoms)) deallocate(atoms)
    end subroutine gr_hutch_mc


    !Updates  gr with new ones.
    !Should be called when hrmc accepts the new position.
    subroutine accept_gr(m, used_data_sets)
        implicit none
        type(model), intent(inout) :: m
        logical, dimension(4), intent(in) :: used_data_sets
        integer i, j, k
        if(used_data_sets(1)) then
            do i=1, mbin_e
                gr_e_sim_cur(i) = gr_e_sim_new(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_e_part_cur(j,k,i) = gr_e_part_new(j,k,i)
                    enddo
                enddo
            enddo
        endif
        if(used_data_sets(2)) then
            do i=1, mbin_n
                gr_n_sim_cur(i) = gr_n_sim_new(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_n_part_cur(j,k,i) = gr_n_part_new(j,k,i)
                    enddo
                enddo
            enddo
        endif
        if(used_data_sets(3)) then
            do i=1, mbin_x
                gr_x_sim_cur(i) = gr_x_sim_new(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_x_part_cur(j,k,i) = gr_x_part_new(j,k,i)
                    enddo
                enddo
            enddo
        endif
    end subroutine accept_gr


    !Retruns gr to the old one.
    !Should be called when hrmc rejects the new position.
    subroutine reject_gr(m,used_data_sets)
        implicit none
        type(model), intent(inout) :: m
        logical, dimension(4), intent(in) :: used_data_sets
        integer i,j,k

        if(used_data_sets(1)) then
            do i=1, mbin_e
                gr_e_sim_new(i) = gr_e_sim_cur(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_e_part_new(j,k,i) = gr_e_part_cur(j,k,i)
                    enddo
                enddo
            enddo
        endif
        if(used_data_sets(2)) then
            do i=1, mbin_n
                gr_n_sim_new(i) = gr_n_sim_cur(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_n_part_new(j,k,i) = gr_n_part_cur(j,k,i)
                    enddo
                enddo
            enddo
        endif
        if(used_data_sets(3)) then
            do i=1, mbin_x
                gr_x_sim_new(i) = gr_x_sim_cur(i)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_x_part_new(j,k,i) = gr_x_part_cur(j,k,i)
                    enddo
                enddo
            enddo
        endif
    end subroutine reject_gr


    !Does the same thing as the subroutine gr_no_hutch, but this subroutine uses hutch algorithm.
    !Written for testing purpose here.
    !Calculates gr with hutch. 
    ! The function is currently not used.
    subroutine gr_hutch(m,used_data_sets,istat)   
        type(model), intent(in) :: m
        integer, intent(inout) :: istat
        !integer status
        integer,  dimension(:), pointer :: atoms
        integer:: nlist        
        integer i, j, k, ig
        real vb, nid, rho
        real xi, yi, zi, xj, yj, zj, xij, yij, zij, R, R2
        real rmax_e,rmax_n,rmax_x
        logical, dimension(4), intent(in) :: used_data_sets
        integer sum_h, sum_no_h

        !******************************************
        rho = m%natoms/(m%lx*m%ly*m%lz)

        do i=1, m%natoms
            sum_h = 0 !added by Feng Yi on 03/10/2009 for debug
            sum_no_h = 0
            !***************************** 
            xi = m%xx%ind(i)
            yi = m%yy%ind(i)
            zi = m%zz%ind(i)

            call hutch_list_3D(m, xi, yi, zi, radius, atoms, istat, nlist)

            do k=1,nlist-1
                if(atoms(k).ne.i)then

                    xj = m%xx%ind(atoms(k))
                    yj = m%yy%ind(atoms(k))
                    zj = m%zz%ind(atoms(k))
            
                    xij = xi - xj
                    yij = yi - yj
                    zij = zi - zj
            
                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))            
                
                    R2 = xij**2+yij**2+zij**2
                    R=sqrt(R2)

                    if(used_data_sets(1)) then
                        rmax_e = del_r_e*mbin_e
                        if(R.lt.rmax_e)then
                          sum_h = sum_h + 1 !debug
                        !IF (i .EQ. 26) THEN
                        !  WRITE(*,*) 'hutch:', atoms(k)
                        !ENDIF !debug

                            ig = int(R/del_r_e) +1
                            gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)+1
                            gr_e_part_new(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_e_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)       
                        endif
                    endif

                    if(used_data_sets(2)) then
                        rmax_n = del_r_n*mbin_n
                        if(R.lt.rmax_n)then
                            ig = int(R/del_r_n) +1
                            gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)+1   
                            gr_n_part_new(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_n_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)       
                        endif
                    endif

                    if(used_data_sets(3)) then
                        rmax_x = del_r_x*mbin_x
                        if(R.lt.rmax_x)then
                            ig = int(R/del_r_x) +1
                            gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)+1
                            gr_x_part_new(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig) = gr_x_part_cur(m%znum_r%ind(i),m%znum_r%ind(atoms(k)),ig)       
                        endif
                    endif
                endif
            enddo !do k
            !WRITE(*,*) 'Size of atoms is: ', size(atoms, 1)
            if(associated(atoms)) DEALLOCATE(atoms) !added by Feng Yi on 03/03/2009

            !*******The following is for debug
            do j=m%natoms + 1, m%natoms
                if(i.ne.j)then
                    !WRITE(*,*) 'j is: ', j !debug
                    xj = m%xx%ind(j)
                    yj = m%yy%ind(j)
                    zj = m%zz%ind(j)
                
                    xij = xi - xj
                    yij = yi - yj
                    zij = zi - zj

                    xij = xij-m%lx*anint(xij/(m%lx))                
                    yij = yij-m%ly*anint(yij/(m%ly))            
                    zij = zij-m%lz*anint(zij/(m%lz))            
        
                    R2 = xij**2+yij**2+zij**2
                    R=sqrt(R2)

                    if(used_data_sets(1)) then
                        rmax_e = del_r_e*mbin_e
                        if(R.lt.rmax_e)  then
                          !if (i .EQ. 26) THEN
                        !   write(*,*) 'no hutch', j
                         ! endif
                          sum_no_h = sum_no_h + 1
                        endif !if R .lt. rmax_e
                    endif !if used_date_sets(1)
                endif !i .ne. j
            enddo !j
            
            !write(*,*) i, 'hutch', sum_h, 'no hutch', sum_no_h,  'difference', sum_h-sum_no_h !debug
            
        enddo !do i

        if(used_data_sets(1)) then
            do i=1,mbin_e
                R = del_r_e*(i-1.0)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_e_sim_cur(i) = gr_e_sim_cur(i) + be(j,k)*gr_e_part_cur(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=912,file="test_tot_e_gr_hutch.txt",form='formatted',status='unknown')
            do i=1, mbin_e
                R = del_r_e*(i)-del_r_e
                vb =((i+1)**3-i**3)*del_r_e**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_e_sim_cur(i) = (gr_e_sim_cur(i)/(m%natoms)-nid)/(R*del_r_e)                
                !write(912,*)R, gr_e_sim_cur(i)
            enddo
            !close(912)
        endif

        if(used_data_sets(2)) then
            do i=1,mbin_n
                R = del_r_n*(i-1.0)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_n_sim_cur(i) = gr_n_sim_cur(i) + bn(j,k)*gr_n_part_cur(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=913,file="test_tot_n_gr_hutch.txt",form='formatted',status='unknown')
            do i=1, mbin_n
                R = del_r_n*(i)-del_r_n
                vb =((i+1)**3-i**3)*del_r_n**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_n_sim_cur(i) = (gr_n_sim_cur(i)/(m%natoms)-nid)/(R*del_r_n)               
                !write(913,*)R, gr_n_sim_cur(i)
            enddo
            !close(913)
        endif

        if(used_data_sets(3)) then
            do i=1,mbin_x
                R = del_r_x*(i-1.0)
                do j=1, m%nelements
                    do k=1, m%nelements
                        gr_x_sim_cur(i) = gr_x_sim_cur(i) + bx(j,k)*gr_x_part_cur(j,k,i) 
                    enddo
                enddo
            enddo

            !open(unit=914,file="test_tot_x_gr_hutch.txt",form='formatted',status='unknown')
            do i=1, mbin_x
                R = del_r_x*(i)-del_r_x
                vb =((i+1)**3-i**3)*del_r_x**3
                nid =(4.0/3.0)*pi*vb*rho
                gr_x_sim_cur(i) = (gr_x_sim_cur(i)/(m%natoms)-nid)/(R*del_r_x)                
                !write(914,*)R, gr_x_sim_cur(i)
            enddo
            !close(914)
        endif

    end subroutine gr_hutch

end module gr_mod
