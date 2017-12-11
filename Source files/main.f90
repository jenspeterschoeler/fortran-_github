program main

    !! The "main" routine of the program  
    !!   
    !! This is where execution starts when running the program
    
    use processor
    use fea

    implicit none

     ! Read model data
     !print*,"Entering Input"
     call input

   	 ! Initialize problem
     !print*,"Entering Initial"
     call initial
     
     ! Calculate displacements
	 select case (antype)
     case('static')
    	call stopwatch('star')
    	call displ
     case('modal')
    	write(*,*) "Input number of Eigen frequencies for solution:"
    	read(*,*)neig
		call stopwatch('star')
    	call eigen_freq
     case('thermal')
     	call thermal
	 case('thermstat')
     
     	call stopwatch('star')
     	call thermal
        call displ
     end select
    ! Close plot window(s)
    call plot( done )
    
end program main
