! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 21.02.2014
! version: 0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use math_funs
    use pre_process_data
    implicit none
    character(len=20)                             :: file_name
    type(tetraElement), dimension(:), allocatable :: tetras
    integer                                       :: npart,i, id
    real(dp), dimension(:,:), allocatable         :: vertices
    real                                          :: time1, time2 ! runtime variables
    
    
    ! user input
    file_name = "ex10.msh"
    
    npart = 8
    call read_mesh_data(file_name, npart, tetras, vertices)
!     deallocate(tetras)
    
!     npart = 2
!     call read_mesh_data(file_name, npart, tetras, vertices)
!     deallocate(tetras)
!      
!     npart = 4
!     call read_mesh_data(file_name, npart, tetras, vertices)
!     deallocate(tetras)
!     
!     npart = 8
!     call read_mesh_data(file_name, npart, tetras, vertices)
!     deallocate(tetras)
!     
!     npart = 16
!     call read_mesh_data(file_name, npart, tetras, vertices)
!     deallocate(tetras) 
     
    
end program raytracing 
 
 
