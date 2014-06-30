! user input of material properties for raytracing and other user input

module rt_properties

	implicit none
    integer, parameter               :: dp = selected_real_kind(p=14) ! precision
    real(dp), parameter              :: pi = 3.14159265358979
    real(dp), parameter              :: sbk = 5.67e-8             ! stefan-boltzmann constant
    character(len=7), parameter      :: objFolder = "../obj/"     ! object folder
    character(len=8), parameter      :: dataFolder = "../data/"   ! data folder
    character(len=11), parameter     :: resFolder = "../results/"   ! results folder
    real(dp)                         :: Etotal = 0.0
    real(dp)                         :: Eleft = 0.0
    integer, parameter               :: npart = 20          ! number of partitions for data input
    integer, parameter               :: nrays = 100000  ! number of rays to be emitted in total    
    integer, parameter               :: nRayPaths = 100  ! number of ray path written out in file (maximal 10000)
    character(len = 100), dimension(2), parameter :: ignoredSurfaces = (/"Iface_Domain1", "Iface_Domain2"/)
    
    ! led setup values
!     character(len=20), parameter     :: data_fname = "led" ! prefix for input data file
!     character(len=20), parameter     :: RT_setup = 'led'  ! select which setup is considered, values are 'tomo' or 'led'
    real(dp), parameter              :: Eemitted = 100.0   ! power to be emitted
    real(dp)                         :: Eblue = 0.0
    real(dp)                         :: Eyellow = 0.0
    real(dp), parameter              :: kappa = 10    ! absorption coefficient
    real(dp), parameter              :: sigma = 10    ! scattering coefficient
    real(dp), dimension(2),parameter :: refracIndices = (/1.9, 1.5/) ! refraction indices
    real(dp), parameter              :: qe = 0.95     ! quantum efficiency (percentag of rays whi get absorbed)
    real(dp), parameter              :: absorptionPercentage = 0.05 ! amount of total absorbtion                                                                                      
!     character(len=100), dimension(1), parameter :: emSurfNames = (/"Iface_Domain2"/)
    
!     ! tomo setup value
    character(len=20), parameter                :: data_fname = "Cube" ! prefix for input data file
!     character(len=20), parameter                :: data_fname = "cube2d" ! prefix for input data file
    character(len=100), dimension(3), parameter :: emSurfNames = (/"zLow_Domain2", "zHigh_Domain2", "Iface_Domain1"/)
!     character(len=100), dimension(2), parameter :: emSurfNames = (/"zLow_Domain2", "zHigh_Domain2"/)
    character(len=20), parameter                :: RT_setup = 'tomo'  ! select which setup is considered, values are 'tomo' or 'led'
    character(len=20), parameter                :: reflection = 'diffuse'  ! select which type of reflection is used (specular or diffuse) for interface
    character(len=30), parameter                :: input = 'faceEmissions_tomo.bin'
    real(dp), parameter                         :: alpha = 0.15    ! hemispherical absorptivity
	real(dp), dimension(2), parameter           :: tbounds = (/373.15, 283.15/) ! temperature for zlow and zhigh walls
    
    ! additional customization in rt_parameters and 
    ! for the 'tomo' setup in rt_tracing subroutine DomainChange
    
end module rt_properties
