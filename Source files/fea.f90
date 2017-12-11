module fea

    !! This module contains procedures that are common to FEM-analysis
    
    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, &
    eigen_freq, mmul, thermal, buildkt, enforcekt, buildthermload, TempElem

contains

!
!----------------------------------------------------------------
!
    subroutine initial
        
        !! This subroutine is mainly used to allocate vectors 
        !! and matrices
        
        use fedata
        use link1
        use plane42rect
        use plane42
		
		! Allocating parameters for half bandwidth
        integer :: e, i, nen, edofMin, edofMax
		integer, dimension(mdim) :: edof
		integer, dimension(4) :: edof2
  
               
! Hint for continuum elements:
!        integer, parameter :: mdim = 8 
!        integer, dimension(mdim) :: edof

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn
		bw = 0
        bw_t = 0
        if (.not. banded) then
            select case (antype)
            case('static')
				allocate (kmat(neqn, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3), & 
                strainthermMat(ne, 3))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
        		c = 0
                strainthermMat = 0
            case('modal')
				allocate (kmat(neqn, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
        		c = 0
            case('thermal')
				allocate (kt(nn, nn))
                allocate (Rt(nn), T(nn), pt(neqn), T_elem(ne))
            case('thermstat')
				allocate (kt(neqn, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn), pt(neqn))
                allocate (kmat(neqn, neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3), &
                strainthermMat(ne, 3), T_elem(ne))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
        		c = 0
              	strainthermMat = 0
            end select 
        else
            ! Find bw
                ! Find coordinates and degrees of freedom
            do e = 1,ne
                edof = 0
                nen = element(e)%numnode
                ! Extrema's are zeroed pr. element
                edofMin = neqn
                edofMax = 0
                ! EDOF is calculated
                do i = 1, nen
                    edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                    edof(2*i)   = 2 * element(e)%ix(i)
                 end do
                 edofMax = MAXVAL(edof)
                 edofMin = MINVAL(edof)
                 ! If the bandwidth exceeds previous 
                 ! elements bw the bw is updated
                 if ((edofMax-edofMin +1) > bw) then 
                     bw = edofMax -edofMin +1
                 end if
            end do

            ! Find bw thermal
			do e = 1,ne
                edof2 = 0
                nen = element(e)%numnode
                ! Extrema's are zeroed pr. element
                edofMin = nn
                edofMax = 0
                ! EDOF is calculated
                do i = 1, nen
                    edof2(i) = element(e)%ix(i)  
                 end do
                 edofMax = MAXVAL(edof2)
                 edofMin = MINVAL(edof2)
                 ! If the bandwidth exceeds previous 
                 ! elements bw the bw is updated
                 if ((edofMax-edofMin +1) > bw_t) then 
                     bw_t = edofMax -edofMin +1
                 end if
            end do
            
            select case (antype)
            case('static')
				allocate (kmat(bw, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3), &
                strainthermMat(ne, 3))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
        		c = 0
                strainthermMat = 0
            case('modal')
				allocate (kmat(bw, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
        		c = 0
            case('thermal')
				allocate (kt(bw_t, nn))
             	allocate (Rt(nn), T(nn), pt(neqn), T_elem(ne))
                pt = 0
            case('thermstat')
				allocate (kt(bw_t, nn))
                allocate (Rt(nn), T(nn))
                allocate (kmat(bw, neqn))
                allocate (p(neqn), d(neqn), p_applied(neqn),&
                 pt(neqn))
                allocate (strain(ne, 3), stress(ne, 3), &
        		stress_vm(ne), principal_stresses(ne, 3), &
                strainthermMat(ne, 3), T_elem(ne))
                ! Initialise parameters
                strain = 0
        		stress = 0
        		stress_vm = 0
        		max_vm = 0
        		principal_stresses = 0
                strainthermMat = 0
                pt = 0
        		c = 0
            end select 
        end if
		
		! intialise parameters found in all antypes
		h2 = 0

 
         
     end subroutine initial
 !
 !-------------------------------------------------------
 !
     subroutine displ
 
         !! This subroutine calculates displacements
 
         use fedata
         use numeth
         use processor
 
         integer :: e, i, j, nen, ix_e(4),max_vm_loc 
         real(wp), dimension(:), allocatable :: plotval, &
         plotval_therm
         real(wp), dimension(mdim) :: xvm_max 
		 ! Build load-vector
         call buildload

         ! Update load vector with thermal nodal loads
         if (antype == 'thermstat') then
           p(1:neqn) = p(1:neqn) + pt(1:neqn)
         end if
         
		 p_applied(1:neqn) = p(1:neqn)
         !print*,"p_applied = ", p_applied
         ! Build stiffness matrix
         call buildstiff
         
		 print*, minval(t), maxval(t), mprop(1)%cond
         ! Remove rigid body modes
         call enforce
         ! print*,"banded = ",banded
         if (banded) then
             call bfactor(kmat)
             call bsolve(kmat, p)
         else
             ! Factor stiffness matrix
             call factor(kmat)
             ! Solve for displacement vector
             call solve(kmat, p)
         end if
         ! Transfer results
         d(1:neqn) = p(1:neqn)
         !print*,'d = ',d
         c = dot_product(d,p_applied)

         ! Recover stress
         call recover
         max_vm = maxval(stress_vm)
		 do i = 1, ne
    			if (stress_vm(i) .eq. max_vm) then
        		max_vm_loc = i
        		exit
    		endif
		 end do
         
         do j = 1, 4
            xvm_max(2*j-1) = x(element(max_vm_loc)%ix(j),1)
            xvm_max(2*j  ) = x(element(max_vm_loc)%ix(j),2)
            !edof(2*j-1) = 2 * element(e_interest)%ix(j) - 1  
            !edof(2*j)   = 2 * element(e_interest)%ix(j)
         end do
			
		 !print*, 'output',mprop(1)%alpha, mprop(2)%alpha,  d(9156)
		 !print*,bw_t,bw,neqn,nn
         
		 call stopwatch("stop") 
     
         ! Output results
         call output
 		 if (antype == 'thermstat') then
          	! Evaluate element temperatures 		 
			call plot(  elements,device=-3, title="Element temperature distribution", &
            legend=.true., eval=T_elem ) ! ,device = -3
         end if
         
         ! Plot deformed shape
         !call plot( deformed+undeformed,device=matlab, &
         ! title = 'Deformed structure' ) ! device=-3
		 !call plot( deformed+undeformed, title = &
         !'Merry Bi-metalic Christmas' )
         
         ! Plot element values
         allocate (plotval(ne))
         do e = 1, ne
             if (element(e)%id == 1) then
                 plotval(e) = stress(e,1)
             else if (element(e)%id == 2) then
                 !! Plots Von mises stress
                 plotval(e) = stress_vm(e)					
				 !! Plots first principal stress
                 !plotval(e) = principal_stresses(e,1)		
                 !! Plots second principal stress
                 !plotval(e) = principal_stresses(e,2)		
             end if
         end do
         !call plot(  elements, device = -3, title="Von mises Stress", &
         ! legend=.true., eval=plotval ) !device = -3
         !call plot ( vectors, p1=principal_stresses(1:ne,1),&
         ! p2=principal_stresses(1:ne,2), &
         !ang=principal_stresses(1:ne,3),title="Principal Stresses") 
         ! print*, "results = ",h2, c, max_vm
     end subroutine displ
 !
 !---------------------------------------------------------------------
 !
     subroutine buildload
 
         !! This subroutine builds the global load vector
 
         use fedata
         use plane42
 
         integer :: i
 ! Hint for continuum elements:
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(mdim) :: re
 
 ! Our own
         integer :: e_interest
         integer :: j
         integer :: nen
         integer :: eface
         real(wp) :: fe
         real(wp) :: thk
 
 ! Initialise hints
         e_interest = 0
         j = 0
         nen = 0
         eface = 0
         fe  = 0
         thk = 0
         edof = 0
         xe = 0
         re = 0
 
         ! Build load vector
         p(1:neqn) = 0
         do i = 1, np            
             select case(int(loads(i, 1)))
             case( 1 )
                 !! Build nodal load contribution
                     !! i,2 is node and i,3 is d.o.f (x or y) 
                     !! and i,4 is the actual load
                 p(2*int(loads(i, 2))-2+int(loads(i,3)))=loads(i,4)             
             case( 2 )
                 !! Build uniformly distributed surface (pressure) 
                 !! load contribution
                 e_interest = int(loads(i,2))
                 eface = int(loads(i,3))
                 fe = loads(i,4)
                 thk  = mprop(element(int(loads(i,2)))%mat)%thk
                 nen = element(e_interest)%numnode
                 do j = 1, nen
                     xe(2*j-1) = x(element(e_interest)%ix(j),1)
                     xe(2*j  ) = x(element(e_interest)%ix(j),2)
                     edof(2*j-1) = 2 * element(e_interest)%ix(j)-1  
                     edof(2*j)   = 2 * element(e_interest)%ix(j)
                 end do
                 call plane42_re(xe,eface,fe,thk,re)
                 !OBS: mdim is hardcoded in fedata
                 do j = 1, mdim
                     p(edof(j)) = p(edof(j)) + re(j)
                 end do

 				
             case default
                 write(*,*)'ERROR in fea/buildload'
                 write(*,*) 'Load type not known'
                 stop
             end select
         end do

     end subroutine buildload
 !
 !----------------------------------------------------------------------
 !
     subroutine buildstiff
 
         !! This subroutine builds the global stiffness matrix from
         !! the local element stiffness matrices
 
         use fedata
         use link1
         use plane42
 
         integer :: e, i, j
         integer :: nen

 ! Hint for system matrix in band form:
 !        integer :: irow, icol
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(mdim, mdim) :: ke 
 ! Hint for modal analysis:
 !        real(wp), dimension(mdim, mdim) :: me 
         real(wp) :: young, area 
 ! Hint for modal analysis and continuum elements:
         real(wp) :: nu, thk		 

         ! Reset stiffness matrix
         kmat = 0
         do e = 1, ne
           
             !! Find coordinates and degrees of freedom
             nen = element(e)%numnode
             do i = 1, nen
                  xe(2*i-1) = x(element(e)%ix(i),1)
                  xe(2*i  ) = x(element(e)%ix(i),2)
                  edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                  edof(2*i)   = 2 * element(e)%ix(i)
             end do
 
             !! Gather material properties and find 
             !! element stiffness matrix
             select case( element(e)%id )
             case( 1 )
                  young = mprop(element(e)%mat)%young
                  area  = mprop(element(e)%mat)%area
                  call link1_ke(xe, young, area, ke)
             case( 2 )
                  young = mprop(element(e)%mat)%young
                  nu  = mprop(element(e)%mat)%nu
                  thk  = mprop(element(e)%mat)%thk
                  call plane42_ke(xe, young, nu, thk, ke)
             end select
 
             ! Assemble into global matrix
             if (.not. banded) then
                 do i = 1, 2*nen
                     do j = 1, 2*nen
                         kmat(edof(i), edof(j)) = kmat(edof(i), &
                         edof(j)) + ke(i, j)
                     end do
                 end do
 ! Hint: Can you eliminate the loops above 
 ! by using a different Fortran array syntax?
             else
             	 do i = 1, 2*nen
                     do j = 1, 2*nen                      
                     	 if (edof(i) >= edof(j)) then
                         	 kmat(edof(i)-edof(j)+1, edof(j))&
                              = kmat(edof(i)-edof(j)+1,edof(j)) &
                              +ke(i, j)
                         end if
                     end do
                 end do
             end if 
         end do
     end subroutine buildstiff
 !
 !-----------------------------------------------------------
 !
     subroutine enforce
 
         !! This subroutine enforces the support boundary conditions
 
         use fedata
 
         integer :: i, j, idof
         real(wp) :: penal
         ! Correct for supports
         if (.not. banded) then
             if (.not. penalty) then
                 do i = 1, nb
                     idof = int(2*(bound(i,1)-1) + bound(i,2))
                     p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) &
                     * bound(i, 3)
                     p(idof) = bound(i, 3)
                     kmat(1:neqn, idof) = 0
                     kmat(idof, 1:neqn) = 0
                     kmat(idof, idof) = 1
                 end do
             else
                 penal = penalty_fac*maxval(kmat)
                 do i = 1, nb
                     idof = int(2*(bound(i,1)-1) + bound(i,2))
                     kmat(idof, idof) = kmat(idof, idof) + penal
                     p(idof) = penal * bound(i, 3)  
                 end do  
             end if
         else
             do i = 1, nb
                 idof = int(2*(bound(i,1)-1) + bound(i,2))
                 !print*,"idof = ",idof
                 !p(1:neqn)=p(1:neqn)-kmat(1:neqn, idof)*bound(i, 3) 
                 !!  doesnt work since kmat is banded
                 p(idof) = 0   ! bound(i, 3)
                 !! column zeroed
                 kmat(1:bw,idof) = 0 
                 !! 1 in the originals diagonal
				 kmat(1,idof) = 1	

                 if(min(bw, idof)>1) then
                   do j = 2,min(bw,idof)
                     kmat(j,idof-j+1) = 0
                   end do
                 end if
				 !do j = 1,bw-1
                ! 	counter = counter + 1
                 !   if (idof - counter <= 0) then
                 !   	exit
                 !   else
                 !   	kmat(j+1,idof-counter) = 0
                 !   end if
                 !end do
             end do
         end if
     end subroutine enforce
 !
 !------------------------------------------------------------------
 !
     subroutine recover
 
         !! This subroutine recovers the element stress, element strain, 
         !! and nodal reaction forces
 
         use fedata
         use link1
         use plane42
 
         integer :: e, i, nen
         integer :: edof(mdim)
         real(wp), dimension(mdim) :: xe, de
         real(wp), dimension(mdim, mdim) :: ke
         real(wp) :: young, area 
 ! Hint for continuum elements:
         real(wp):: nu, sigma_vm, he
         real(wp), dimension(3) :: estrain, estress,&
         principal_stresses_e, estraintherm

         ! Reset force vector
         p = 0
 
         do e = 1, ne
           
             ! Find coordinates etc...
             nen = element(e)%numnode
             do i = 1,nen
                 xe(2*i-1) = x(element(e)%ix(i), 1)
                 xe(2*i)   = x(element(e)%ix(i), 2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
                 de(2*i-1) = d(edof(2*i-1))
                 de(2*i)   = d(edof(2*i))
             end do
 
             ! Find stress and strain
             select case( element(e)%id )
             case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
                 p(edof(1:2*nen)) = p(edof(1:2*nen)) &
                 + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                 call link1_ss(xe, de, young, estress, estrain)
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain
             case( 2 )
             	 young = mprop(element(e)%mat)%young
				 nu = mprop(element(e)%mat)%nu
                 estraintherm = strainthermMat(e,1:3)
                 call plane42_ss(xe, de, young, nu, estraintherm, &
                 estress, estrain, principal_stresses_e, sigma_vm, he)
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain
                 stress_vm(e) = sigma_vm
                 h2 = h2 + he
                 principal_stresses(e, 1:3) = principal_stresses_e          
             end select
         end do
         h2 = h2/ne
     end subroutine recover
 !
 !------------------------------------------------------------------
 !
     subroutine mmul(invector, outvector)
		 !! This subroutine elementwise multiplication of mass matrix
         !! with eigenvector
         use fedata
         use plane42
         
 		 real(wp), dimension(:), intent(in) :: invector
         real(wp), dimension(:), intent(out) :: outvector
         
         integer :: e, i, j, ng
         integer :: nen

 ! Hint for system matrix in band form:
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
 ! Hint for modal analysis:
         real(wp), dimension(mdim, mdim) :: me 
 ! Hint for modal analysis and continuum elements:
         real(wp) :: dens, thk, eta, xi, detjac
         real(wp) :: w(2), xi_vec(2), eta_vec(2), n(2,8),&
         bmat(3,8), jac(2,2), out_e(8), in_e(8)
 		 
		 ng = 2

         outvector = 0
         do e = 1, ne
           	 me = 0
             in_e = 0
             out_e = 0
             ! Find coordinates and degrees of freedom
             nen = element(e)%numnode
             do i = 1, nen
                  xe(2*i-1) = x(element(e)%ix(i),1)
                  xe(2*i  ) = x(element(e)%ix(i),2)
                  edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                  edof(2*i)   = 2 * element(e)%ix(i)
             end do
             ! Elementwise invector
 			 do i = 1,2*nen
               in_e(i) = invector(edof(i))
             end do
             
             ! Gather material properties and find element mass matrix
             thk  = mprop(element(e)%mat)%thk
             dens  = mprop(element(e)%mat)%dens
             eta_vec(1) = -1/sqrt(3.0)
             eta_vec(2) = 1/sqrt(3.0)
             xi_vec(1) = -1/sqrt(3.0)
             xi_vec(2) = 1/sqrt(3.0)
           	 w(1) = 1
             w(2) = 1  
             do i = 1,ng
                   xi = xi_vec(i)
                   do j = 1,ng
                       eta = eta_vec(j)
                       call plane42_shape(xe, xi, eta,&
                        n, bmat, jac, detjac)
                       me = me + dens*w(i)*w(j)&
                       *matmul(transpose(n),n)*detjac*thk
                   end do
             end do
             out_e = matmul(me,in_e)
             ! Global outvector
             do i = 1,2*nen
               outvector(edof(i)) = outvector(edof(i)) + out_e(i)
             end do
	 	 end do
     end subroutine mmul

 !
 !---------------------------------------------------------------
 !
     subroutine eigen_freq
 
         !! This subroutine calculates eigenvalues
 
         use fedata
         use numeth
         use processor

         !! Mass matrix and eigenvalue calc
         real(wp) :: xp(neqn), yp(neqn), yp_old(neqn), xp_old(neqn), &
         yp_old_temp(neqn), Dre(neqn, neig), lambda, rp, omega_n(neig),&
         Z(neqn, neig), cj(neig), Z_in(neqn), Dre_in(neqn), Z_out(neqn)
         integer :: p_incr, i, pmax, idof, j, q
         real(wp), dimension(:), allocatable :: plotval
		 
		 
         yp = 0
         xp_old = 0
         yp_old = 0
		 p = 0

         ! Build stiffness matrix
         call buildstiff

         ! Remove rigid body modes
         call enforce
         if (.not. banded) then
             ! Factor stiffness matrix
             call factor(kmat)
         else
             call bfactor(kmat)
         end if
         
   		 Dre = 0
         do i = 1,neig
            ! Guess for initial eigenvector
            xp = 0.001
            Z = 0
            ! Build initial Yp
            call mmul(xp,yp)
            ! Build Z_j
            do j = 1, i-1
              	Dre_in = Dre(:,j)
				call mmul(Dre_in,Z_out)
                Z(:,j) = Z_out
            end do
            !Inverse iteration loop
            pmax = 1000
            do p_incr = 1,pmax
               xp_old = xp
               yp_old = yp
               ! Enforce BC's on yp
               do q = 1, nb
                   idof = int(2*(bound(q,1)-1) + bound(q,2))
                   yp_old(idof) = bound(q, 3)
               end do
               yp_old_temp = yp_old
               call bsolve(kmat,yp_old)
               xp = yp_old
               yp_old = yp_old_temp
               ! Compute cj
			   do j = 1, i-1
               		Z_in = Z(:,j)
               		cj(j) = dot_product(xp,Z_in)
                    xp = xp - cj(j)*Dre(:,j)
               end do              
               call mmul(xp,yp)
               rp = sqrt(dot_product(xp,yp))
               yp = yp/rp
               if (sqrt(dot_product((xp-xp_old),(xp-xp_old)))&
                 /sqrt(dot_product(xp,xp)) < 1e-12 ) then
                   exit
               end if
            end do
            Dre(:,i) = xp/rp         
            lambda = dot_product(xp,yp_old)/rp**2
         	omega_n(i) = sqrt(lambda)
         end do
         call stopwatch("stop")      
		 print*,'omega_n = ',omega_n
         do j = 1,neig
         	D = Dre(:,j)
		 	call plot( eigenmode, device=-3,&
             eigenfreq=omega_n(j)&
            , eigenvector=D, title="Eigen plot" , & 
			tend=5.0_wp, tdelta=100.0_wp )
         end do
     end subroutine eigen_freq 
     
 !
 !---------------------------------------------------------------
 !
     
     subroutine thermal
 
         !! This subroutine calculates eigenvalues
 
         use fedata
         use numeth
         use processor
         use plane42
         use plane41		          

		 real(wp), dimension(:), allocatable :: plotval_therm
         integer :: e, i, nen, ix_e(4)
         real(wp) :: Tmax, T_elemMax
         call buildthermload

         call buildkt
         !do i = 1,nn
         !	 print*,'line',i
         !   print*,kt(i,1:nn)
         !end do
		 call enforcekt
         
         if (banded) then
             call bfactor(kt)
             call bsolve(kt, Rt)
         else
             ! Factor stiffness matrix
             call factor(kt)
             ! Solve for displacement vector
             call solve(kt, Rt)
         end if
         ! Transfer results
         T(1:nn) = Rt(1:nn)
         Tmax = maxval(T)
		 
		 ! Calculate element temperatures
		 call TempElem
         
		 T_elemMax = maxval(T_elem)


         ! Output to terminal
         !print*,'T = ', T
         !print*,'T_nodal_max = ', Tmax
         !print*,'T_elem_max = ',T_elemMax
		 if (antype == 'thermstat') then
		 	call recoverTherm
         end if	 

		 if ( antype == 'thermal') then
            ! Plot element temperatures         
         	call plot(  elements, title="Element temperature distribution", &
          	legend=.true., eval=T_elem)
         end if
     end subroutine thermal

 !
 !----------------------------------------------------------------------
 !
     subroutine buildkt
 
         !! This subroutine builds the global stiffness matrix from
         !! the local element stiffness matrices
 
         use fedata
         use plane42
         use plane41
 
         integer :: e, i, j
         integer :: nen, eface

 ! Hint for system matrix in band form:
 !        integer :: irow, icol
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(4, 4) :: ce, he
         real(wp) :: cond, conv, thk	 

         ! Reset stiffness matrix
         if (.not. banded) then
             kt = 0
         else
             kt = 0
         end if

         do e = 1, ne           
             !! Find coordinates and degrees of freedom
             nen = element(e)%numnode
             do i = 1, nen
                  xe(2*i-1) = x(element(e)%ix(i),1)
                  xe(2*i  ) = x(element(e)%ix(i),2)
                  edof(i) = element(e)%ix(i)
 !                 edof(2*i)   = 2 * element(e)%ix(i)
             end do
 
             !! Gather material properties and find 
             !! thermal stiffness matrix
             select case( element(e)%id )
             case( 1 )
                  write(*,*)'Case: Truss Elements - Not defined for thermal analysis'
             case( 2 )
                  he = 0
                  ce = 0
                  cond = mprop(element(e)%mat)%cond
                  thk  = mprop(element(e)%mat)%thk
                  call plane41_ce(xe, cond, thk, ce)
                  ! Convection from fluid flow
                  do i = 1,ncp                   
                  	if (int(convloads(i,2)) == e ) then
                    	conv = convloads(i, 4)
                        eface = convloads(i, 3)
                  		call plane41_he(xe, eface, conv, thk, he)
                  	end if
                  end do
             end select
 
             ! Assemble into global matrix
             if (.not. banded) then
                 do i = 1, nen
                     do j = 1, nen                       	 
                         kt(edof(i), edof(j)) = kt(edof(i),&
                          edof(j)) + ce(i, j) + he(i, j)                     
                     end do
                 end do               

             else
             	 do i = 1, nen
                     do j = 1, nen                      
                     	 if (edof(i) >= edof(j)) then
                         	 kt(edof(i)-edof(j)+1, edof(j)) = &
                             kt(edof(i)-edof(j)+1, edof(j))&
                              + ce(i,j)&
                              + he(i, j)
                         end if
                     end do
                 end do
             end if 
         end do

     end subroutine buildkt

 !
 !-----------------------------------------------------------
 !
     subroutine enforcekt
 
         !! This subroutine enforces the boundary conditions
 
         use fedata
 
         integer :: i, idof!, j, counter!is for 0 & 1's method (banded)
         real(wp) :: penal!, bwcolumn(nn)!is for 0 & 1's method (banded)
         ! Correct for supports
         if (.not. banded) then
             if (.not. penalty) then
                 do i = 1, ntb
					 idof = int(tbound(i,1))
                     Rt(1:nn) = Rt(1:nn)-kt(1:nn, idof)*tbound(i,3)  
                     Rt(idof) = tbound(i,3) 
                     kt(1:nn, idof) = 0
                     kt(idof, 1:nn) = 0
                     kt(idof, idof) = 1.0
                 end do
             else
                 penal = 1e9*maxval(kt)
                 do i = 1, ntb
                     idof = int(tbound(i,1))
                     kt(idof, idof) = kt(idof, idof) + penal
                     Rt(idof) = Rt(idof) + penal * tbound(i,3)
                 end do  
             end if
         else
         	if (.not. penalty) then
!$$$$$$              do i = 1, ntb
!$$$$$$                  idof = int(tbound(i,1))
!$$$$$$                  bwcolumn(1:nn) = 0
!$$$$$$                  bwcolumn(idof:int(idof+bw_t-1)) = kt(1:bw_t, idof)
!$$$$$$                  print*,'nydsa'
!$$$$$$                  print*,bwcolumn
!$$$$$$                  !if(min(bw, idof)>1) then
!$$$$$$                  !   do j = 2,min(bw,idof)
!$$$$$$                  !       print*,'j = ',j,'   ', 'idof = ',idof
!$$$$$$                  !       bwcolumn(idof-j+1) = kt(j,idof-j+1)
!$$$$$$                  !        print*,bwcolumn(idof-j+1)
!$$$$$$                  !       end do
!$$$$$$                  !end if
!$$$$$$                  !p(1:neqn)=p(1:neqn)-kmat(1:neqn, idof)*bound(i, 3) 
!$$$$$$                  !!  doesnt work since kmat is banded
!$$$$$$                  Rt(1:nn) = Rt(1:nn) - bwcolumn(1:nn)*tbound(i,3)
!$$$$$$                  !Rt(idof) = tbound(i, 3)
!$$$$$$                  counter = 0
!$$$$$$                  !! column zeroed
!$$$$$$                  kt(1:bw_t,idof) = 0 
!$$$$$$                  !! 1 in the originals diagonal
!$$$$$$                  kt(1,idof) = 1 
!$$$$$$ 
!$$$$$$                  if(min(bw_t, idof)>1) then
!$$$$$$                    do j = 2,min(bw_t,idof)
!$$$$$$                      kt(j,idof-j+1) = 0
!$$$$$$                    end do
!$$$$$$                  end if
!$$$$$$                  !do j = 1,bw-1
!$$$$$$                  !  counter = counter + 1
!$$$$$$                  !   if (idof - counter <= 0) then
!$$$$$$                  !      exit
!$$$$$$                  !   else
!$$$$$$                  !      kmat(j+1,idof-counter) = 0
!$$$$$$                  !   end if
!$$$$$$                  !end do
!$$$$$$             end do
             else
                 penal = 1e9*maxval(kt)
                 do i = 1, ntb
                     idof = int(tbound(i,1))
                     kt(1, idof) = kt(1, idof) + penal
                     Rt(idof) = Rt(idof) + penal * tbound(i,3)
                 end do  
             end if
         end if
     end subroutine enforcekt
 !
 !------------------------------------------------------------------
 !

     subroutine buildthermload
 
         !! This subroutine builds the global load vector
 
         use fedata
         use plane42
         use plane41
 
         integer :: i, e, j, z
 ! Hint for continuum elements:
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(4) :: rhe, rbe, rqe
 
 ! Our own
         integer :: e_interest
         integer :: nen
         integer :: eface
         real(wp) :: Tfl, conv, fb, Quni
         real(wp) :: thk

 ! Initialise hints
         e_interest = 0
         j = 0
         nen = 0
         eface = 0
         Tfl  = 0
         conv = 0
         thk = 0
         edof = 0
         xe = 0
         rhe = 0
         rbe = 0
         fb = 0


         ! Build thermal load vector
         Rt(1:nn) = 0
		 do i = 1, ncp
                 !! Build load from convection
                 e_interest = int(convloads(i,2))
                 eface = int(convloads(i,3))
                 conv = convloads(i,4)
                 Tfl = convloads(i,5)
                 thk  = mprop(element(e_interest)%mat)%thk
                 nen = element(e_interest)%numnode
                 do j = 1, nen
                     xe(2*j-1) = x(element(e_interest)%ix(j),1)
                     xe(2*j  ) = x(element(e_interest)%ix(j),2)
                     edof(j) = element(e_interest)%ix(j)
                 end do
                 call plane41_rhe(xe, eface, conv, thk, Tfl, rhe)
                 do j = 1, 4
                     Rt(edof(j)) = Rt(edof(j)) + rhe(j)
                 end do
          end do

          do i = 1, nhp
                 !! Build load from heat flux
                 e_interest = int(fluxloads(i,2))
                 eface = int(fluxloads(i,3))
                 fb = fluxloads(i,4)
                 thk  = mprop(element(e_interest)%mat)%thk
                 nen = element(e_interest)%numnode
                 do j = 1, nen
                     xe(2*j-1) = x(element(e_interest)%ix(j),1)
                     xe(2*j  ) = x(element(e_interest)%ix(j),2)
                     edof(j) = element(e_interest)%ix(j)
                 end do
                 call plane41_rbe(xe, eface, thk, fb, rbe)
                 do j = 1, 4
                     Rt(edof(j)) = Rt(edof(j)) + rbe(j)
                 end do
          end do
          
          if (.not. nq == 0) then
          	do i = 1,nq
            	Quni = qgen(i,2)
                ! All heat generations are put on all elements
          		do e = 1, ne
				 	thk  = mprop(element(e)%mat)%thk
                 	nen = element(e)%numnode
                 	do j = 1, nen
                     	xe(2*j-1) = x(element(e)%ix(j),1)
                     	xe(2*j  ) = x(element(e)%ix(j),2)
                     	edof(j) = element(e)%ix(j)
                 	end do
                 	call plane41_rqe(xe, thk, Quni, rqe)
                 	do j = 1, 4
                    	 Rt(edof(j)) = Rt(edof(j)) + rqe(j)
                 	end do
            	end do
             end do
          end if
          if (.not. nqe == 0) then
          	do z = 1,nqe
            	 Quni = qegen(z,2)
                 ! All heat generations are put on all elements
				 thk  = mprop(element(int(qegen(z,1)))%mat)%thk
                 nen = element(int(qegen(z,1)))%numnode
                 do j = 1, nen
                     xe(2*j-1) = x(element(int(qegen(z,1)))%ix(j),1)
                     xe(2*j  ) = x(element(int(qegen(z,1)))%ix(j),2)
                     edof(j) = element(int(qegen(z,1)))%ix(j)
                 end do
                 call plane41_rqe(xe, thk, Quni, rqe)
                 do j = 1, 4
                   	 Rt(edof(j)) = Rt(edof(j)) + rqe(j)
                 end do

             end do
          end if
         
     end subroutine buildthermload

 !
 !------------------------------------------------------------------
 !
     subroutine recoverTherm
 
         !! This subroutine recovers the element stress, element strain, 
         !! and nodal reaction forces
 
         use fedata
         use plane42
         use plane41
 
         integer :: e, i, nen, j, m, k, z
         integer :: edof(mdim)
         real(wp), dimension(mdim) :: xe
         real(wp) :: young, nu, thk, alpha
         real(wp) :: cmat(3,3)
         real(wp) :: straintherm(3), pt_e(8), Te_e
         real(wp) :: detjac, bmat(3,8), jac(2,2), n(2,8), xi, eta, &
         gauss(2), wg
		 integer :: ix_e(4), ng
         

         ! gaussian points
         gauss(1) 	= (-1)/((3.0)**0.5) !xi - node1
         gauss(2) 	= (1)/((3.0)**0.5) !eta- node1
         wg = 1 ! weight function
         ng = 2
         
         ! Reset thermal force vector
         pt = 0
         do e = 1, ne 
           	 pt_e = 0     
             ! Find coordinates etc...
             nen = element(e)%numnode

             do z = 1,nen
                 xe(2*z-1) = x(element(e)%ix(z), 1)
                 xe(2*z)   = x(element(e)%ix(z), 2)
                 edof(2*z-1) = 2 * element(e)%ix(z) - 1
                 edof(2*z)   = 2 * element(e)%ix(z)
				 ix_e(z) = element(e)%ix(z) 
             end do
             ! Find Thermal strain
            
			 Te_e = T_elem(e)
             young = mprop(element(e)%mat)%young
			 nu = mprop(element(e)%mat)%nu
             thk = mprop(element(e)%mat)%thk
             alpha = mprop(element(e)%mat)%alpha
             call plane41_ssTherm(Te_e, young, nu, alpha, cmat, straintherm)
             do i = 1,ng
                xi = gauss(i)
                do j = 1,ng
                   eta = gauss(j)
                   call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
                   pt_e = pt_e+wg*wg*matmul(matmul(transpose(bmat),&
                   cmat),straintherm)*detjac
                end do
             end do
             do k = 1,8
             	pt(edof(k)) = pt(edof(k)) + pt_e(k)
             end do
             strainthermMat(e, 1:3) = straintherm
         end do
     end subroutine recoverTherm

 !
 !------------------------------------------------------------------
 !
     subroutine TempElem
 	 	 use fedata
         use plane41
		
         integer :: e,i, edof(4), nen
         real(wp), dimension(mdim) :: xe
         real(wp) :: T_e
         real(wp) :: T_nodal(4)

         do e = 1,ne
            nen = element(e)%numnode
            do i = 1,nen
                 xe(2*i-1) = x(element(e)%ix(i), 1)
                 xe(2*i)   = x(element(e)%ix(i), 2)
                 edof(i) = element(e)%ix(i)
            end do
            T_nodal = T(edof)
			call plane41_Te(xe,T_nodal, T_e)
            T_elem(e) = T_e
         end do

     end subroutine TempElem
 !
 !---------------------------------------------------------------------
 !
 
 end module fea
 

