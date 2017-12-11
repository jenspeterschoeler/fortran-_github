module plane41

    !! This module contains subroutines specific to the plane41 element 
    !!        
    !! The plane41 element has 4 nodes. Each node 
    !! has 1 degree-of-freedom, namely, temperature
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
    public ::  plane41_ce, plane41_shapeTherm, &
    plane41_he, plane41_rhe, plane41_rbe, &
    plane41_rqe, plane41_ssTherm, plane41_Te
   

contains

 subroutine plane41_ce(xe, cond, thk, ce)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: cond
            !! Conductivity for this element
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
        real(wp), dimension(:,:), intent(inout) :: ce
            !! Conductivity matrix
        real(wp) :: bmatTherm(2,4), nTherm(4), jacTherm(2,2), Imat(2,2) 
        real(wp) :: eta, xi, detjacTherm
		integer :: i,j, ng
		real(wp) :: w(2), xi_vec(2), eta_vec(2)        

        !Build unit matrix
		Imat = 0
        Imat(1,1) = cond
        Imat(2,2) = cond
        
		ng = 2
        ! print*,'ce not zeroed in plane_41_ce'
		! ce = 0
          
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
                    call plane41_shapeTherm(xe, xi,&
                     eta, nTherm, bmatTherm, jacTherm,&
                      detjacTherm)
                    ce = ce + w(i)*w(j)*matmul(matmul&
                    (transpose(bmatTherm),Imat),&
                    bmatTherm)*detjacTherm*thk                               
                end do
            end do
        else
          	print*,'error in plane42_ce'
          	print*,'you need to add code for gauss points =/= 2'
		end if

    end subroutine plane41_ce

!
!-----------------------------------------------------------
!    

subroutine plane41_he(xe, eface, conv, thk, he)
    
        !! This subroutine computes the boundary convection matrix due
        !! to convection heat flux (flux is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where convection happens
        real(wp), intent(in) :: conv
            !! Convection heat transfer coefficient
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])
        real(wp), dimension(:,:), intent(inout) :: he
            !! Boundary convection matrix
            
        real(wp) :: gauss(2), xi, &
        eta,bmatTherm(2,4),detjacTherm, &
        jacTherm(2,2),nTherm(4),nmul(1,4),nmult(4,1)
        real(wp) :: L_12, L_23, L_34, L_41
        ! Define Gauss points
        integer :: wg, i,ng
        
		! gaussian points
        gauss(1) 	= (1)/((3.0)**0.5) !xi - node1
        gauss(2) 	= (-1)/((3.0)**0.5) !eta- node1
        wg = 1 ! weight function
        ng = 2
	
		! he = 0 not zeroed here in the function but outisde in the loop because 
        ! 		thermal stiffness has to be higher when a corner node is considered 
		if (eface == 1) then
        	L_12 = sqrt((xe(3)-xe(1))**2+(xe(4)-xe(2))**2)
        	eta = (-1)
			do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi, eta,&
                 nTherm, bmatTherm, jacTherm, detjacTherm)
                nmul(1,1:4) = nTherm(1:4)
                nmult(1:4,1) = nTherm(1:4)
          		he =he+ wg * thk * conv &
                * matmul(nmult,nmul)*L_12*0.5
        	end do
            
        elseif (eface == 2) then
        	L_23 = sqrt((xe(5)-xe(3))**2+(xe(6)-xe(4))**2)
			xi = 1
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, nTherm, &
                bmatTherm, jacTherm, detjacTherm)
                nmul(1,1:4) = nTherm(1:4)
                nmult(1:4,1) = nTherm(1:4)
          		he =he+ wg * thk * conv &
                * matmul(nmult,nmul)*L_23*0.5
        	end do

        elseif (eface == 3) then
        	L_34 = sqrt((xe(7)-xe(5))**2+(xe(8)-xe(6))**2)
			eta=1
        	do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, nTherm,&
                 bmatTherm, jacTherm, detjacTherm)
                nmul(1,1:4) = nTherm(1:4)
                nmult(1:4,1) = nTherm(1:4)
          		he =he+ wg * thk * conv &
                * matmul(nmult,nmul)*L_34*0.5
        	end do

        elseif (eface == 4) then
        	L_41 = sqrt((xe(7)-xe(1))**2+(xe(8)-xe(2))**2)
        	xi=(-1)
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, nTherm,&
                 bmatTherm, jacTherm, detjacTherm)
                nmul(1,1:4) = nTherm(1:4)
                nmult(1:4,1) = nTherm(1:4)
          		he =he+ wg * thk * conv &
                * matmul(nmult,nmul)*L_41*0.5
        	end do

        endif
    
    end subroutine plane41_he


!
!-----------------------------------------------------------
!    

subroutine plane41_rhe(xe, eface, conv, thk, Tfl, rhe)
    
        !! This subroutine computes the boundary convection vector
        
        integer, intent(in) :: eface
            !! Element face where convection happens
        real(wp), intent(in) :: conv
            !! Convection heat transfer coefficient
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), intent(in) :: Tfl
            !! Bulk temperature of convecting fluid
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])
        real(wp), dimension(:), intent(out) :: rhe
            !! Boundary convection vector
        real(wp) :: L_12, L_23, L_34, L_41
            
        real(wp) :: gauss(2), xi, &
        eta,bmatTherm(2,4),detjacTherm, jacTherm(2,2), nTherm(4)
        ! Define Gauss points
        integer :: wg, i,ng
        
		! gaussian points
        gauss(1) 	= (1)/((3.0)**0.5) !xi - node1
        gauss(2) 	= (-1)/((3.0)**0.5) !eta- node1
        wg = 1 ! weight function
        ng = 2
	
		rhe=0
		if (eface == 1) then
            L_12 = sqrt((xe(3)-xe(1))**2+(xe(4)-xe(2))**2)
        	eta = (-1)
			do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, &
                nTherm, bmatTherm, jacTherm, detjacTherm)
          		rhe =rhe+ wg * thk * conv * Tfl &
                * nTherm*L_12*0.5
        	end do
            
        elseif (eface == 2) then
        	L_23 = sqrt((xe(5)-xe(3))**2+(xe(6)-xe(4))**2)
			xi = 1
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, &
                nTherm, bmatTherm, jacTherm, detjacTherm)
          		rhe =rhe+ wg * thk * conv * Tfl &
                * nTherm*L_23*0.5
        	end do

        elseif (eface == 3) then
            L_34 = sqrt((xe(7)-xe(5))**2+(xe(8)-xe(6))**2)
			eta=1
        	do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, &
                nTherm, bmatTherm, jacTherm, detjacTherm)
          		rhe =rhe+ wg * thk * conv * Tfl &
                * nTherm*L_34*0.5
        	end do

        elseif (eface == 4) then
            L_41 = sqrt((xe(7)-xe(1))**2+(xe(8)-xe(2))**2)
        	xi=(-1)
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi, eta, &
                nTherm, bmatTherm, jacTherm, detjacTherm)
          		rhe =rhe+ wg * thk * conv * Tfl &
                * nTherm*L_41*0.5
        	end do

        endif
    
    end subroutine plane41_rhe

!
!-----------------------------------------------------------
!    

subroutine plane41_rqe(xe, thk, Quni, rqe)
    
        !! This subroutine computes the heat generation nodal load vector
        
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), intent(in) :: Quni
            !! Uniform heat generation of body
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])
        real(wp), dimension(:), intent(out) :: rqe
            !! Boundary convection vector
            
        real(wp) :: xi_vec(2),eta_vec(2) ,xi, &
        eta,bmatTherm(2,4),detjacTherm, jacTherm(2,2)
		real(wp), dimension(4) :: nTherm
        ! Define Gauss points
        integer :: w(2), i, j,ng
        
        ng = 2
		rqe = 0
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
                    call plane41_shapeTherm(xe, xi,&
                     eta, nTherm, bmatTherm, &
                     jacTherm, detjacTherm)
					rqe = rqe + w(i)*w(j)&
                    *Quni*thk*detjacTherm*nTherm
                end do
            end do
        else
          	print*,'error in plane41_rqe'
          	print*,'you need to add code for gauss points =/= 2'
		end if
    
    end subroutine plane41_rqe

!
!----------------------------------------------------------------
!
    subroutine plane41_shapeTherm(xe, xi, eta, &
    nTherm, bmatTherm, jacTherm, detjacTherm)

        !! This subroutine computes the shape 
        !! functions, strain displacement matrix,
        !! the Jacobian matrix and the determinant 
        !! of the Jacobian matrix for a specific Gauss point

        real(wp), intent(in) :: xi
            !! Natural coordinate xi
        real(wp), intent(in) :: eta
            !! Natural coordinate eta
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed
            !! configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(out) :: nTherm
            !! Shape function matrix
        real(wp), dimension(:,:), intent(out) :: bmatTherm
        	!! Strain displacement matrix
        real(wp), dimension(:,:), intent(out) :: jacTherm
        	!! Jacobian matrix
        real(wp), intent(out) :: detjacTherm
        	!! Determinant of jacobian matrix "shape factor"

        real(wp) :: NbarTherm(2, 8), gamma(2,2)
        real(wp) :: Nbar_fact, jac_fact, gamma_fact, n_fact
		
		
		! Setup NbarTherm matrix
		Nbar_fact = 0.25
        NbarTherm = 0.0      
        ! Unassigned elements are all zero
		NbarTherm(1,1) = (-1.+eta)*Nbar_fact !N1,xi
        NbarTherm(1,2) = (1.-eta)*Nbar_fact  !N2,xi
        NbarTherm(1,3) = (1.+eta)*Nbar_fact  !N3,xi
        NbarTherm(1,4) = (-1.-eta)*Nbar_fact !N4,xi
        NbarTherm(2,1) = (-1.+xi)*Nbar_fact  !N1,eta
        NbarTherm(2,2) = (-1.-xi)*Nbar_fact  !N2,eta
        NbarTherm(2,3) = (1.+xi)*Nbar_fact   !N3,eta
        NbarTherm(2,4) = (1.-xi)*Nbar_fact   !N4,eta       
        
		! Setup Jacobian matrix
        jac_fact = 0.25
        jacTherm = 0.0
        
        jacTherm(1,1) = ((-1+eta)*xe(1)+(1-eta)*xe(3)&
        +(1+eta)*xe(5)+(-1-eta)*xe(7))*jac_fact
        jacTherm(2,1) = ((-1+xi)*xe(1)+(-1-xi)*xe(3)&
        +(1+xi)*xe(5)+(1-xi)*xe(7))*jac_fact
        jacTherm(1,2) = ((-1+eta)*xe(2)+(1-eta)*xe(4)&
        +(1+eta)*xe(6)+(-1-eta)*xe(8))*jac_fact
        jacTherm(2,2) = ((-1+xi)*xe(2)+(-1-xi)*xe(4)&
        +(1+xi)*xe(6)+(1-xi)*xe(8))*jac_fact


        ! Determinant of Jacobian
		detjacTherm = jacTherm(1,1)*jacTherm(2,2)&
        -jacTherm(2,1)*jacTherm(1,2)

        ! Setup gamma (inverted Jacobian)
        gamma_fact = 1.0/detjacTherm
        gamma = 0.0
        gamma(1,1) = jacTherm(2,2)*gamma_fact
        gamma(2,1) = -jacTherm(2,1)*gamma_fact
        gamma(1,2) = -jacTherm(1,2)*gamma_fact
        gamma(2,2) = jacTherm(1,1)*gamma_fact


		! Compute strain-displacement matrix
        bmatTherm = matmul(gamma,NbarTherm)  

        ! Setup nTherm
        n_fact = 0.25
        nTherm = 0.0
        nTherm(1) = (1-xi)*(1-eta)*n_fact
        nTherm(2) = (1+xi)*(1-eta)*n_fact
        nTherm(3) = (1+xi)*(1+eta)*n_fact
        nTherm(4) = (1-xi)*(1+eta)*n_fact
        
    end subroutine plane41_shapeTherm

!
!-----------------------------------------------------------
!    

subroutine plane41_rbe(xe, eface, thk, fb, rbe)
    
        !! This subroutine computes the thermal load due
        !! to convection.

        integer, intent(in) :: eface
            !! Element face where convection happens
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), intent(in) :: fb
            !! Heat flux
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])
        real(wp), dimension(:), intent(out) :: rbe
            !! Heat flux vector
        real(wp) :: L_12, L_23, L_34, L_41
            
        real(wp) :: gauss(2), xi, &
        eta,bmatTherm(2,4),detjacTherm, jacTherm(2,2), nTherm(4)

        ! Define Gauss points
        integer :: wg, i,ng
        
		! gaussian points
        gauss(1) 	= (1)/((3.0)**0.5) !xi - node1
        gauss(2) 	= (-1)/((3.0)**0.5) !eta- node1
        wg = 1 ! weight function
        ng = 2

		rbe=0
		if (eface == 1) then
            L_12 = sqrt((xe(3)-xe(1))**2+(xe(4)-xe(2))**2)
        	eta = (-1)
			do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi, &
                eta, nTherm, bmatTherm, jacTherm, detjacTherm)
          		rbe =rbe+ wg * thk * fb &
                * nTherm*L_12*0.5
        	end do
            
        elseif (eface == 2) then
        	L_23 = sqrt((xe(5)-xe(3))**2+(xe(6)-xe(4))**2)
			xi = 1
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi, eta,&
                 nTherm, bmatTherm, jacTherm, detjacTherm)
          		rbe =rbe+ wg * thk * fb &
                * nTherm*L_23*0.5
        	end do

        elseif (eface == 3) then
            L_34 = sqrt((xe(7)-xe(5))**2+(xe(8)-xe(6))**2)
			eta=1
        	do i = 1,ng
            	xi = gauss(i)
                call plane41_shapeTherm(xe, xi,&
                 eta, nTherm, bmatTherm, jacTherm, detjacTherm)
          		rbe =rbe+ wg * thk * fb &
                * nTherm*L_34*0.5
        	end do

        elseif (eface == 4) then
            L_41 = sqrt((xe(7)-xe(1))**2+(xe(8)-xe(2))**2)
        	xi=(-1)
        	do i = 1,ng
            	eta = gauss(i)
                call plane41_shapeTherm(xe, xi&
                , eta, nTherm, bmatTherm, jacTherm, detjacTherm)
          		rbe =rbe+ wg * thk * fb &
                * nTherm*L_41*0.5
        	end do

        endif
    end subroutine plane41_rbe

	!
	!-----------------------------------------------------------
	!   
    
    subroutine plane41_ssTherm(Te_e, young,&
     nu, alpha, cmat, straintherm)

        !! This subrotuine computes nodal forces from
        !! thermal expansion

        real(wp), intent(in) :: Te_e
            !! Element average temperature 
        real(wp), intent(in) :: alpha
            !! Thermal expansion coefficient
        real(wp), intent(in) :: nu
            !! Poisson's ratio
        real(wp), intent(in) :: young
            !! youngs modulus
        real(wp), dimension(:), intent(out) :: straintherm
        	!! Nodal forces from thermal expansion

        real(wp), dimension(:,:), intent(out) :: cmat
        real(wp) :: fact
	
        
        ! Build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1.-nu**2.)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1.-nu)/2.

        straintherm(1) = alpha*Te_e
        straintherm(2) = straintherm(1)
        straintherm(3) = 0


    end subroutine plane41_ssTherm
!
!-----------------------------------------------------------
!    

subroutine plane41_Te(xe, T_nodal, T_e)
    
        !! This subroutine computes the heat generation nodal load vector
        

            !! Thickness of this element

            !! Uniform heat generation of body
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed 
            !! configuration (see also [[plane42_ke]])

        real(wp), dimension(:), intent(in) :: T_nodal
        	!! Nodal temperature
        real(wp),  intent(out) :: T_e
            !! Boundary convection vector
        real(wp) :: gauss(2), xi, &
        eta,bmatTherm(2,4),detjacTherm, jacTherm(2,2), nTherm(4)
           
		
		xi = 0
        eta = 0

        call plane41_shapeTherm(xe, xi, eta,&
         nTherm, bmatTherm, jacTherm, detjacTherm)
		T_e = dot_product(nTherm,T_nodal)
    end subroutine plane41_Te

    
end module plane41


