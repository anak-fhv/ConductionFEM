! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 09.05.2014
! version: 0.6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use pre_process
    use tracing

    implicit none
    integer  :: i
    real(dp) :: dummy
    
    ! user input
    ! all user input is set in rt_properties and used as global variables
    
    ! intialize random generator
    dummy = myRandom(2803) 
    
    call start_preprocessing
!     write(*,*) tetraData(36348)%vertexIDs
!     write(*,*) tetraData(36348)%neighbors(:,1)
!     write(*,*) tetraData(36348)%neighbors(:,2)
!     write(*,*)
!     write(*,*) tetraData(36220)%vertexIDs
!     write(*,*) tetraData(36220)%neighbors(:,1)
!     write(*,*) tetraData(36220)%neighbors(:,2)
    call start_tracing
    
    ! just some output
    open(unit=101, file=objFolder//"absorbed.res", action='write',status='new')   
	write(101,'(e14.6)') (absorbed(i), i =1, size(absorbed))
    close(unit=101)
!     
!     write(*,*) "fraction of power absorbed:", sum(absorbed)/Etotal
    
end program raytracing 