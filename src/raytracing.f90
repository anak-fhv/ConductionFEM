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
    type(rayContainer)                               :: ray
    type(tetraElement), dimension(:), allocatable    :: tetraData
    type(emissionSurface), dimension(:), allocatable :: emSurf
    real(dp), dimension(:,:), allocatable            :: vertices
    integer                                          :: npart, k
    real(dp)                                         :: dummy, t1, t2
    
    ! user input
    file_name = "sphere.dat"
    npart = 12
    emSurfNames = ['xLow'] ! right now this requires to state the correct dimension in above declaration       
    
    ! some internal bookkeeping
    dummy = myRandom(2803)
    
    ! perform pre-processing
    ! depending on the suffix of the input file either the pre-processing routine ...
    if (index(file_name,".msh") > 0) then
        call read_mesh_data(file_name, emSurfNames, npart, tetraData, vertices, emSurf)
        call WriteData(tetraData, vertices, emSurf,file_name)
    else
        ! ... or an exisiting file is read
        call ReadData(tetraData, vertices, emSurf, file_name)
    end if
    
    call cpu_time(t1)
    ! create ray
    do k = 1,10
        call CreateRay(tetraData, vertices, emSurf(1), ray)
        call TraceRay(tetraData, vertices, ray)
    end do
    call cpu_time(t2)
    
    write(*,*) "run time", t2-t1
   
end program raytracing 