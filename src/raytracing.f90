! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 08.04.2014
! version: 0.4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use pre_process_data
    use tracing

    implicit none
    character(len=100)                               :: file_name
    character(len=100), dimension(1)                 :: emSurfNames
    type(tetraElement), dimension(:), allocatable    :: tetraData
    type(emissionSurface), dimension(:), allocatable :: emSurf
    real(dp), dimension(:,:), allocatable            :: vertices
    integer                                          :: npart, nrays
    
    ! user input
    file_name = "sphere"
    npart = 4
    emSurfNames = ['xLow'] ! right now this requires to state the correct dimension in above declaration       
    nrays = 100000
    
    call start_preprocessing(file_name, emSurfNames, npart, tetraData, vertices, emSurf)
    call start_tracing(tetraData, vertices, emSurf, file_name, nrays)
   
end program raytracing 