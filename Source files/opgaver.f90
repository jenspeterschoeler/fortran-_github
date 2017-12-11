module opgaver

    !! This module contains procedures that are specifically used to solve assignments
    
    implicit none
    save
    private
    public :: opg63

contains

!
!--------------------------------------------------------------------------------------------------
!

    subroutine opg63
    	! This subroutine is used to find the elements of interest in question 6.3
		use fedata
		real(wp) :: B_r_cond
		integer :: e, i, A_node, B_node, B_node_right
        
		! Compliance v. size plot data
		 c = dot_product(d,p_applied)
		 B_r_cond = 1+esize_width(1)
         B_node_right = 0
         print*,"esize_width = ",esize_width
         print*,"esize_height = ",esize_height
		 print*,"B_r_cond = ",B_r_cond
         do i = 1,nn
           	!print*,"x = ",x(i,1),"y = ",x(i,2)
            if (x(i,1) == 0 .AND. x(i,2) == 1) then
                A_node = i
                print*,"A_node = ", A_node
            elseif (x(i,1) == 1 .AND. x(i,2) == 1) then
                B_node = i
                print*,"B_node = ", B_node
            elseif (x(i,1) == B_r_cond .AND. x(i,2) == 1) then
            	B_node_right = i
                print*,"B_node_right = ", B_node_right
            end if
         end do
         	
         do e = 1,ne
            do i = 1,4
              	if(element(e)%ix(i) == A_node) then
                	element_interestA = e
                    print*,"element_interestA = ",element_interestA
                end if
            end do
         end do    

		 do e = 1,ne
            if (element(e)%ix(1) == B_node .OR. element(e)%ix(2) == B_node & 
               .OR. element(e)%ix(3) == B_node .OR. element(e)%ix(4) == B_node) then
				do i = 1,4
                	if(element(e)%ix(i) == B_node_right) then
                    	element_interestB = e
                        print*,'element_interestB = ',element_interestB
                    end if
                end do
         	end if
         end do        

    end subroutine opg63
 
 end module opgaver