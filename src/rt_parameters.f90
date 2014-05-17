! RayTracing modules defining types and global (internal) variables
! also contains functions which require user input or change

module rt_parameters

	use rt_properties
    implicit none
    
    type :: tetraElement
        integer, dimension(4)    :: vertexIds=0   ! ID ofvertex vertices of the tetra element
        
        integer, dimension(4,2)  :: neighbors=0   ! information about of neighboring elements
        ! In above array the row index relates to the 4 faces of the element. 
        ! The 1st column contains the ID of the neighboring element or surface. In case of
        ! a surface the ID = ID_of_surface + Number_of_Tetraelements.
        ! The 2nd column contains the face of the neighboring element which is identical to
        ! the face given by the row index. In case a surface is the neighbor the entry is -1.
        
        !real(dp), dimension(4,4) :: shape_funcs ! shape functions
        
        integer                  :: domain=0      ! to which domain the tetra belongs
    end type tetraElement
    
    type :: emissionSurface
        character(len=100)                   :: name ! name of emission surface
        real(dp), dimension(:), allocatable  :: area ! cumsum of area of the faces on the surface of emission
        real(dp)                             :: totalarea ! totalarea of surface
        integer, dimension(:,:), allocatable :: elemData ! (:,1) number of tetra-element      
                                                         ! (:,2) face of tetra-element which is on the surface
    end type
    
    type :: rayContainer  ! contains information on the traced ray
        real(dp), dimension(3) :: point, direction     ! current point and its direction
        integer                :: tetraID              ! current tetraeder index
        integer                :: faceID               ! current face index
        real(dp)               :: length = 0.0_dp      ! distance travelled
        real(dp)               :: power = 1.0_dp       ! current power of ray
        real(dp)               :: wavelength = 1.0_dp  ! wavelength
    end type
    
    ! some gloabl variables
    real(dp), dimension(:), allocatable              :: absorbed  ! field containing info about absorptio
    real(dp), dimension(:,:), allocatable            :: vertices  ! filed of all vertices
    type(tetraElement), dimension(:), allocatable    :: tetraData ! type for tetraeder information
    type(emissionSurface), dimension(:), allocatable :: emSurf    ! emission surfaces
    
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
    
    
    ! select an emission surface
     integer function emsIDfun()
	    
	    real(dp)             :: r, totalarea
	    integer              :: i
	    
	    ! the code below selects an emission surface based on
	    ! the overall area of all emission surfaces
	    ! get totalarea of all emission surfaces
	    totalarea = sum(emSurf%totalarea)
	    r = myRandom(0)
		
		do emsIDfun = 1,size(emSurf)
			if (r <= sum(emSurf(1:emsIDfun)%totalarea)/totalarea) exit
	    end do
	    		 
	end function emsIDfun
	
	
	! determine power of a single ray
	real(dp) function raypowerfun(temp, area)
		
		real(dp), intent(in) :: temp  ! temperature (in K)
		real(dp), intent(in) :: area  ! area of emission face 
		
		! setup if power is independent of location and temperature
! 		raypowerfun = Etotal/nrays
		
		! setup if temperature-dependence exist
		raypowerfun = eta*sbk*temp**4*area
		
	end function raypowerfun
	
    
end module rt_parameters







        
    
    
    