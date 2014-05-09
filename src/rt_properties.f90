! user input of material properties for raytracing and other user input
! does also include wrapper for random numbers

module rt_properties

	implicit none
    integer, parameter               :: dp = selected_real_kind(p=14) ! precision
    real(dp), parameter              :: pi = 3.14159265358979
    character(len=7), parameter      :: objFolder = "../obj/"   ! object folder
    character(len=8), parameter      :: dataFolder = "../data/" ! data folder
    real(dp), parameter              :: Etotal = 100.0 ! total energy
    real(dp), parameter              :: kappa = 0.3    ! absorption coefficient
    real(dp), parameter              :: sigma = 0.4    ! scattering coefficient
    real(dp), dimension(2),parameter :: refracIndices = (/1.9, 1.5/) ! refraction indices for different domains
                                                                     ! must be of same size as domains! 
     
    contains
    
    ! phase function for scattering
	subroutine PhaseFunction(theta, psi)
	
		real(dp), intent(out) :: theta, psi
		 
		! isotropic case
	    theta = acos(1.0_dp-2.0_dp*myRandom(0))
	    psi = 2.0_dp*pi*myRandom(0)		
        
    end subroutine PhaseFunction
    
    
    ! wrapper for random numbers
    function myRandom(iflag)
        use ifport           ! use intel random numbers           
        
        real(dp) :: myRandom
        integer  :: iflag
        
        if (iflag > 0) call srand(iflag) ! seeds random number
        myRandom = drand(0)
        
    end function myRandom
    
end module rt_properties