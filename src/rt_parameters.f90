! RayTracing modules defining types and global (internal) variables
! also contains functions which require user input or change

module rt_parameters

	use rt_properties
    implicit none
    
    type :: tetraElement
        integer, dimension(4)   :: vertexIds=0   ! ID ofvertex vertices of the tetra element
        
        integer, dimension(4,2) :: neighbors=0   ! information about of neighboring elements
        ! In above array the row index relates to the 4 faces of the element. 
        ! The 1st column contains the ID of the neighboring element or surface. In case of
        ! a surface the ID = ID_of_tetra.
        ! The 2nd column contains the face of the neighboring element which is identical to
        ! the face given by the row index. In case a surface is the neighbor the entry is -ID_of_surface.
        
        !real(dp), dimension(4,4) :: shape_funcs ! shape functions
        
        integer                 :: domain = 0    ! to which domain the tetra belongs
!         real(dp)                :: powerAvail = 0.0  ! power to be emitted from the tetra (used for tomo-setup)
!         integer                 :: nrays = 0     ! number of rays emitted
!         integer                 :: erays = 0     ! number of rays emitted so far
!         real(dp)                :: ppray = 0.0   ! power per ray
    end type tetraElement
    
    type :: emissionSurface
        character(len=100)                   :: name ! name of emission surface
        real(dp), dimension(:), allocatable  :: value ! cumsum of value of the faces on the surface of emission
        real(dp)                             :: totalvalue ! totalvalue of surface
        integer, dimension(:,:), allocatable :: elemData ! (:,1) number of tetra-element      
                                                         ! (:,2) face of tetra-element which is on the surface
        integer                              :: originalID ! id used for tetraeder connections             
        real(dp)                             :: power  ! total power to be emitted from the surface       
        integer, dimension(:,:), allocatable :: rays   ! number of rays emitted from (:,1)and absorbed by (:,2) the face
    end type
    
    type :: rayContainer  ! contains information on the traced ray
        real(dp), dimension(3) :: point, direction     ! current point and its direction
        integer                :: tetraID              ! current tetraeder index
        integer                :: faceID               ! current face index
        real(dp)               :: length = 0.0_dp      ! distance travelled
        real(dp)               :: power = 1.0_dp       ! current power of ray
        real(dp)               :: wavelength = 1.0_dp  ! wavelength
        logical                :: colorchange = .false. ! flag whether color change has happened 
    end type
    
    type :: runstatistic ! statistics for empirical runs
	    real(dp) :: maxvalue, mean, var, rms, logvar,logmean,logmeannormal, logvarnormal, dist
	    integer  :: entries
    end type
    
    ! some gloabl variables
    real(dp), dimension(:), allocatable              :: powerNodal  ! power values distributed to nodes
    real(dp), dimension(:,:), allocatable            :: vertices    ! filed of all vertices
    type(tetraElement), dimension(:), allocatable    :: tetraData   ! type for tetraeder information
    type(emissionSurface), dimension(:), allocatable :: emSurf      ! emission surfaces
    real(dp), dimension(:,:), allocatable              :: source      ! emission for faces on interface for tomo-based setup
    real(dp), dimension(:,:), allocatable            :: spectrumB   ! blue spectrum data for led setup
    real(dp), dimension(:,:), allocatable            :: spectrumY   ! yellow spectrum data for led setup
    
    contains
    ! wrapper for random numbers
    function myRandom(iflag)
        use ifport           ! use intel random numbers           
        
        real(dp) :: myRandom
        integer  :: iflag
        
        if (iflag > 0) call srand(iflag) ! seeds random number
        myRandom = drand(0)
        
    end function myRandom
    
    
    ! phase function for scattering
	subroutine PhaseFunction(theta, psi)
	
		real(dp), intent(out) :: theta, psi
		 
		! isotropic case
	    theta = acos(1.0_dp-2.0_dp*myRandom(0))
	    psi = 2.0_dp*pi*myRandom(0)		
        
    end subroutine PhaseFunction
	
	! get mean free pathlength
	real(dp) function GetPathLength()
	
		! only meaningful for LED setup
	    ! kappa and sigma are globally defined material properties
		GetPathLength = 1.0_dp/(kappa+sigma)*log(1.0_dp/myRandom(0))
	    
	end function GetPathLength
	
    
end module rt_parameters







        
    
    
    