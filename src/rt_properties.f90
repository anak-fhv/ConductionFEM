! user input of material properties for raytracing and other user input

module rt_properties

	implicit none
    integer, parameter               :: dp = selected_real_kind(p=14) ! precision
    real(dp), parameter              :: pi = 3.14159265358979
    real(dp), parameter              :: sbk = 5.67e-8             ! stefan-boltzmann constant
    character(len=7), parameter      :: objFolder = "../obj/"     ! object folder
    character(len=8), parameter      :: dataFolder = "../data/"   ! data folder
    character(len=11), parameter     :: resFolder = "../results/"   ! results folder
    character(len=20), parameter     :: data_fname = "RPC_2d_128" ! prefix for input data file
!     character(len=20), parameter     :: data_fname = "led" ! prefix for input data file
!     real(dp)                         :: Etotal = 100.0 ! total energy
    real(dp)                         :: Etotal = 0.0
    real(dp)                         :: Eleft = 0.0
    real(dp)                         :: Eblue = 0.0
    real(dp)                         :: Eyellow = 0.0
    real(dp), parameter              :: kappa = 20    ! absorption coefficient
    real(dp), parameter              :: sigma = 20    ! scattering coefficient
    real(dp), parameter              :: alpha = 0.15    ! hemispherical absorptivity
    real(dp), dimension(2),parameter :: refracIndices = (/1.9, 1.5/) ! refraction indices for different domains
                                                                     ! must be correspond to domain IDs, last entry is outside of model !
    ! names of surfaces from which rays can be emitted
    ! dimension attribute must be set correctly                                                                  
    character(len=100), dimension(3), parameter      :: emSurfNames = (/"zLow_Domain2", "zHigh_Domain2", "Iface_Domain2"/) 
!     character(len=100), dimension(1), parameter      :: emSurfNames = (/"Iface_Domain1"/)                                                                         
    ! parameter below lists surfaces which do not count as neighbor, but from which rays can be emitted                                                                  
    character(len = 100), dimension(2), parameter    :: ignoredSurfaces = (/"Iface_Domain1", "Iface_Domain2"/)                                                                 
    integer, parameter :: npart = 10          ! number of partitions for data input
    integer, parameter :: nrays = 1000000    ! number of rays to be emitted in total
    character(len=20), parameter :: RT_setup = 'tomo'  ! select which setup is considered, values are 'tomo' or 'led'
    integer, parameter :: nRayPaths = 100  ! number of ray path written out in file (maximal 10000)
    
    ! additional customization in rt_parameters and 
    ! for the 'tomo' setup in rt_tracing subroutine DomainChange
    
end module rt_properties