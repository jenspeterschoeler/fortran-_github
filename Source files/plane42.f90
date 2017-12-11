module plane42

    !! This module contains subroutines specific to the plane42 element 
    !!        
    !! The plane42 element has 4 nodes. Each node 
    !! has 2 degrees-of-freedom, namely, displacement 
    !! along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below.
    !! The figure also shows how the element edges are labelled. 
    !! For example, edge 1 is the edge between node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o 
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!             
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc  
    !! `E1` = element face 1, `E2` = element face 2, etc  
    
    use types
    implicit none
    save
    
    private
    public :: plane42_ke, plane42_re, plane42_ss, plane42_shape!, & 
    !plane42_ce, plane42_shapeTherm, plane42_he, plane42_rhe, plane42_rbe, &
    !plane42_rqe, plane42_ssTherm

contains

    subroutine plane42_ke(xe, young, nu, thk, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in 
            !! undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix
        real(wp) :: cmat(3,3), fact, bmat(3,8), &
        jac(2,2), N(2,8)
        real(wp) :: eta, xi, detjac
		integer :: i,j, ng
		real(wp) :: w(2), xi_vec(2), eta_vec(2)
        
        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2
		
		ng = 2
		ke = 0
          
        if(ng == 2) then
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
                    call plane42_shape(xe, xi, &
                     eta, n, bmat, jac, detjac)
                    ke = ke + w(i)*w(j)*matmul(transpose(bmat)&
                    ,matmul(cmat,bmat))*detjac*thk
                end do
            end do
        else
          	print*,'error in plane42_ke'
          	print*,'you need to add code for gauss points =/= 2'
		end if

    end subroutine plane42_ke
!
!---------------------------------------------------------------
!

    subroutine plane42_ss(xe, de, young, nu, estraintherm, estress,&
     estrain, principal_stresses_e, sigma_vm, he)

        !! This subrotuine computes the element stress and strain 
        !! (((The location inside the element
        !! where stress and and strain is evaluated,
        !! is defined inside the subroutine))).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu 
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed
            !! configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 
            !! along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 
            !! along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(in)  :: estraintherm
        	!! element thermal strain
        real(wp), intent(out) :: sigma_vm
        	!! Initialisation of VM 
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)

        real(wp), dimension(:), intent(out) :: principal_stresses_e
        real(wp), intent(out) :: he
        
        real(wp) :: bmat(3, 8), cmat(3, 3), n(2,8), jac(2,2), detjac, mstrain(3)

        ! Variables for extrapolation of stresses 
        !! from gauss points to nodes
        !real(wp) :: gstressx(4), gstressy(4),&
        ! gstressxy(4), gstress(3), gstrain(3)
        !real(wp) :: gstrainx(4), gstrainy(4),&
        ! gstrainxy(4), nstressx(4), nstressy(4)
        !real(wp) :: nstressxy(4), extra_mat(4,4) 

		! Principal stresses and VM
		real(wp) :: psi, sigma_1, sigma_2, fact

		! hardcoded to the centroid/ >>upper<< of the element
		real(wp) :: xi, eta      
        !real(wp), dimension(4) ::  xi_vec, eta_vec!For gauss point exptrapolation
        !integer :: i!For gauss point exptrapolation
        

!$$$$$$         xi_vec(1) = -1/(sqrt(3.0))
!$$$$$$         xi_vec(2) = 1/(sqrt(3.0))
!$$$$$$         xi_vec(3) = 1/(sqrt(3.0))
!$$$$$$         xi_vec(4) = -1/(sqrt(3.0))
!$$$$$$   
!$$$$$$         eta_vec(1) = -1/(sqrt(3.0))
!$$$$$$         eta_vec(2) = -1/(sqrt(3.0))
!$$$$$$         eta_vec(3) = 1/(sqrt(3.0))
!$$$$$$         eta_vec(4) = 1/(sqrt(3.0))
        
!$$$$$$         ! Build extrapolation matrix for gauss to nodal
!$$$$$$         extra_mat = 0
!$$$$$$         extra_mat(1,1) = 1+0.5*sqrt(3.0)
!$$$$$$         extra_mat(1,2) = -0.5
!$$$$$$         extra_mat(1,3) = 1-0.5*sqrt(3.0)
!$$$$$$         extra_mat(1,4) = -0.5
!$$$$$$         extra_mat(2,1) = -0.5
!$$$$$$         extra_mat(2,2) = 1+0.5*sqrt(3.0)
!$$$$$$         extra_mat(2,3) = -0.5
!$$$$$$         extra_mat(2,4) = 1-0.5*sqrt(3.0)
!$$$$$$         extra_mat(3,1) = 1-0.5*sqrt(3.0)
!$$$$$$         extra_mat(3,2) = -0.5
!$$$$$$         extra_mat(3,3) = 1+0.5*sqrt(3.0)
!$$$$$$         extra_mat(3,4) = -0.5
!$$$$$$         extra_mat(4,1) = -0.5
!$$$$$$         extra_mat(4,2) = 1-0.5*sqrt(3.0)
!$$$$$$         extra_mat(4,3) = -0.5
!$$$$$$         extra_mat(4,4) = 1+0.5*sqrt(3.0)
            

        ! Build constitutive matrix (plane stress)
        cmat = 0.
        fact = young/(1.-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1.-nu)*0.5
        mstrain = 0.
        estrain = 0.
        estress = 0.
        
!$$$$$$         do i = 1, 4   
!$$$$$$            xi = xi_vec(i)
!$$$$$$            eta = eta_vec(i)
!$$$$$$            ! Build strain-displacement matrix
!$$$$$$            call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
!$$$$$$    
!$$$$$$            ! Compute gauss strain
!$$$$$$            gstrain = matmul(bmat, de)
!$$$$$$            gstrainx(i) = gstrain(1)
!$$$$$$            gstrainy(i) = gstrain(2)
!$$$$$$            gstrainxy(i) = gstrain(3) 
           
!$$$$$$            ! Compute gauss stress
!$$$$$$            gstress = matmul(cmat, gstrain)
!$$$$$$            gstressx(i) = gstress(1)
!$$$$$$            gstressy(i) = gstress(2)
!$$$$$$            gstressxy(i) = gstress(3) 
!$$$$$$         end do
        
!$$$$$$         ! Extrapolate from gauss to nodal stress
!$$$$$$         nstressx = matmul(extra_mat,gstressx)
!$$$$$$         nstressy = matmul(extra_mat,gstressy)
!$$$$$$         nstressxy = matmul(extra_mat,gstressxy)
!$$$$$$         
!$$$$$$         ! Element stresses are averaged 
!$$$$$$         estress(1) = sum(nstressx)/4.0
!$$$$$$         estress(2) = sum(nstressy)/4.0
!$$$$$$         estress(3) = sum(nstressxy)/4.0


        xi = 0.
        eta = 0.
        ! Build strain-displacement matrix
        call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
   
        ! Compute element strain
        estrain = matmul(bmat, de)
        mstrain = estrain-estraintherm
        !print*,'mstrain=',mstrain
        estress = matmul(cmat,mstrain)
        
        ! Compute principal stress and direction
        sigma_1 = 0.5*(estress(1)+estress(2))&
        +sqrt((0.5*(estress(1)-estress(2)))**2+estress(3)**2)
        sigma_2 = 0.5*(estress(1)+estress(2))&
        -sqrt((0.5*(estress(1)-estress(2)))**2+estress(3)**2)
        sigma_vm = sqrt(sigma_1**2+sigma_2**2-sigma_1*sigma_2)
        if (sigma_1-sigma_2 == 0.) then
          psi = 0.
        else
        	psi = atan2(-2.*estress(3)/(sigma_1-sigma_2),&
        	(estress(1)-estress(2))/(sigma_1-sigma_2))*0.5
        end if

        principal_stresses_e(1) = sigma_1
        principal_stresses_e(2) = sigma_2
        principal_stresses_e(3) = psi        


        he = 4*detjac
        



    end subroutine plane42_ss

!
!----------------------------------------------------------------
!
    subroutine plane42_shape(xe, xi, eta, n, bmat, jac, detjac)

        !! This subroutine computes the shape functions, strain displacement matrix,
        !! the Jacobian matrix and the determinant of the Jacobian matrix for a specific Gauss point

        real(wp), intent(in) :: xi
            !! Natural coordinate xi
        real(wp), intent(in) :: eta
            !! Natural coordinate eta
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed
            !! configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:,:), intent(out) :: n
            !! Shape function matrix
        real(wp), dimension(:,:), intent(out) :: bmat
        	!! Strain displacement matrix
        real(wp), dimension(:,:), intent(out) :: jac
        	!! Jacobian matrix
        real(wp), intent(out) :: detjac
        	!! Determinant of jacobian matrix "shape factor"

        real(wp) :: L(3,4), Nbar(4, 8), gamma(2,2), gammabar(4,4)
        real(wp) :: Nbar_fact, jac_fact, gamma_fact, n_fact
		
		! Setup L matrix
        L = 0.0
        ! Unassigned elements are all zero
        L(1,1) = 1.0
        L(2,4) = 1.0
        L(3,2) = 1.0
        L(3,3) = 1.0

        ! Setup N
        n_fact = 0.25
        n = 0.0
        n(1,1) = (1-xi)*(1-eta)*n_fact
        n(2,2) = (1-xi)*(1-eta)*n_fact
        n(1,3) = (1+xi)*(1-eta)*n_fact
        n(2,4) = (1+xi)*(1-eta)*n_fact
        n(1,5) = (1+xi)*(1+eta)*n_fact
        n(2,6) = (1+xi)*(1+eta)*n_fact
        n(1,7) = (1-xi)*(1+eta)*n_fact
        n(2,8) = (1-xi)*(1+eta)*n_fact
		
		! Setup Nbar matrix
		Nbar_fact = 0.25
        Nbar = 0.0      
        ! Unassigned elements are all zero
        Nbar(1,1) = (-1+eta)*Nbar_fact  !N1,xi
        Nbar(2,1) = (-1+xi)*Nbar_fact   !N1,eta
        Nbar(3,2) = (-1+eta)*Nbar_fact  !N1,xi
        Nbar(4,2) = (-1+xi)*Nbar_fact   !N1,eta
        Nbar(1,3) = (1-eta)*Nbar_fact   !N2,xi
        Nbar(2,3) = (-1-xi)*Nbar_fact   !N2,eta
        Nbar(3,4) = (1-eta)*Nbar_fact	!N2,xi
        Nbar(4,4) = (-1-xi)*Nbar_fact   !N2,eta
        Nbar(1,5) = (1+eta)*Nbar_fact   !N3,xi
        Nbar(2,5) = (1+xi)*Nbar_fact    !N3,eta
        Nbar(3,6) = (1+eta)*Nbar_fact   !N3,xi
        Nbar(4,6) = (1+xi)*Nbar_fact    !N3,eta
        Nbar(1,7) = (-1-eta)*Nbar_fact  !N4,xi
        Nbar(2,7) = (1-xi)*Nbar_fact    !N4,eta
        Nbar(3,8) = (-1-eta)*Nbar_fact  !N4,xi
        Nbar(4,8) = (1-xi)*Nbar_fact    !N4,eta
        
		! Setup Jacobian matrix
        jac_fact = 0.25
        jac = 0.0
        
        jac(1,1) = ((-1+eta)*xe(1)+(1-eta)*xe(3)&
        +(1+eta)*xe(5)+(-1-eta)*xe(7))*jac_fact
        jac(2,1) = ((-1+xi)*xe(1)+(-1-xi)*xe(3)&
        +(1+xi)*xe(5)+(1-xi)*xe(7))*jac_fact
        jac(1,2) = ((-1+eta)*xe(2)+(1-eta)*xe(4)&
        +(1+eta)*xe(6)+(-1-eta)*xe(8))*jac_fact
        jac(2,2) = ((-1+xi)*xe(2)+(-1-xi)*xe(4)&
        +(1+xi)*xe(6)+(1-xi)*xe(8))*jac_fact


        ! Determinant of Jacobian
		detjac = jac(1,1)*jac(2,2)-jac(2,1)*jac(1,2)
        !print*,'detjac = ',detjac
        
        ! Setup gamma
        gamma_fact = 1.0/detjac
        gamma = 0.0
        gamma(1,1) = jac(2,2)*gamma_fact
        gamma(2,1) = -jac(2,1)*gamma_fact
        gamma(1,2) = -jac(1,2)*gamma_fact
        gamma(2,2) = jac(1,1)*gamma_fact

        ! Setup gammabar
        gammabar = 0.0
        ! Unassigned elements are all zero
        gammabar(1,1) = gamma(1,1)
        gammabar(2,1) = gamma(2,1)
        gammabar(1,2) = gamma(1,2)
        gammabar(2,2) = gamma(2,2)
        gammabar(3,3) = gamma(1,1)
        gammabar(4,3) = gamma(2,1)
        gammabar(3,4) = gamma(1,2)
        gammabar(4,4) = gamma(2,2)
        
		! Compute strain-displacement matrix
        bmat = matmul(matmul(L,gammabar),Nbar)
        
        
    end subroutine plane42_shape
!
!-----------------------------------------------------------
!    

subroutine plane42_re(xe, eface, fe, thk, re)
    

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied
        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at
            !! element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at 
            !! element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: nface(2,1), gauss(2), xi, &
        eta,bmat(3,8),detjac,jac(2,2),n(2,8), re_temp(8,1)
        ! Define Gauss points
        integer :: wg, i,ng
        
		! gaussian points
        gauss(1) 	= (1)/((3.0)**0.5) !xi - node1
        gauss(2) 	= (-1)/((3.0)**0.5) !eta- node1
        wg = 1 ! weight function
        ng = 2
	
		re=0
        re_temp = 0
		if (eface == 1) then
        	eta = (-1)
			do i = 1,ng
            	xi = gauss(i)
                call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
                nface(1,1) = -jac(1,2)
            	nface(2,1) = jac(1,1)
          		re =re+ wg * thk * fe &
                * matmul(transpose(n), nface)
        	end do
            
        elseif (eface == 2) then
			xi = 1
        	do i = 1,ng
            	eta = gauss(i)
                call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
				nface(1,1) = -jac(2,2)
            	nface(2,1) = jac(2,1)
          		re_temp =re_temp+ wg * thk * fe &
                * matmul(transpose(n), nface)
        	end do

        elseif (eface == 3) then
			eta=1
        	do i = 1,ng
            	xi = gauss(i)
                call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
            	nface(1,1) = jac(1,2)
            	nface(2,1) = -jac(1,1)
          		re_temp =re_temp+ wg * thk * fe &
                * matmul(transpose(n), nface)
        	end do

        elseif (eface == 4) then
        xi=(-1)
        do i = 1,ng
            	eta = gauss(i)
                call plane42_shape(xe, xi, eta, n, bmat, jac, detjac)
            	nface(1,1) = jac(2,2)
            	nface(2,1) = -jac(2,1)
          		re_temp =re_temp+ wg * thk * fe &
                * matmul(transpose(n), nface)
        	end do

        endif
    	re = re_temp(1:8,1)
    end subroutine plane42_re
!
!---------------------------------------------------------------
!

    
end module plane42


