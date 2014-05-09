! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 08.04.2014
! version: 0.4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use pre_process
    use tracing

    implicit none
    character(len=100)                               :: file_name
    character(len=100), dimension(1)                 :: emSurfNames    
    type(emissionSurface), dimension(:), allocatable :: emSurf
    integer                                          :: npart, nrays, i
    
    ! user input
    file_name = "sphere"
    npart = 4
    emSurfNames = ['xLow'] ! right now this requires to state the correct dimension in above declaration       
    nrays = 10000
    
    call start_preprocessing(file_name, emSurfNames, npart, emSurf)
    call start_tracing(emSurf, file_name, nrays)
    
    ! just some output
    open(unit=101, file=objFolder//"absorbed.res", action='write',status='new')   
    do i = 1,size(absorbed)
	    if (absorbed(i) /= 0.0_dp) write(101,'(i8,1x,e14.6)') i, absorbed(i)
    end do
    close(unit=101)
    
end program raytracing 