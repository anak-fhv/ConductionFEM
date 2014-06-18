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
!     dummy = myRandom(2803) 
    
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
    open(unit=101, file=resFolder//"absorbed.res", status='replace')   
	write(101,'(e14.6)') (powerNodal(i), i =1, size(powerNodal))
    close(unit=101)
!   
    write(*,*) "power emitted:", Etotal  
    write(*,*) "power absorbed:", sum(powerNodal)
    write(*,*) "power leaving:", Eleft
    write(*,*) "difference:", Etotal + Eleft - sum(powerNodal)
    write(*,*)
    write(*,*) "Eblue:", Eblue
    write(*,*) "Eyellow:", Eyellow
    write(*,*) "max.absorbed", maxval(powerNodal)
    write(*,*) "min.absorbed", minval(powerNodal)
    write(*,*) "mean.absorbed", sum(powerNodal, abs(powerNodal).gt.0)/count(abs(powerNodal).gt.0)
    
    open(unit=102, file=resFolder//"ems.res", status='replace')   
	write(102,'(i8,1x,i8)') (emSurf(3)%rays(i,:), i =1, size(emSurf(3)%rays,1))
    close(unit=102)
    
!     do i=1,size(emSurf(3)%rays,1)
! 	    write(*,*) i, emSurf(3)%rays(i,1)
!     end do
!     
end program raytracing 