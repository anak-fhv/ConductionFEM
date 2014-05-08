! RayTracing modules
! author: Steffen Finck
! contact: steffen.finck@fhv.at


module rt_constants

    implicit none
    integer, parameter          :: dp = selected_real_kind(p=14)
    real(dp), parameter         :: pi = 3.14159265358979
    real(dp), parameter         :: trans = 0.75_dp
    character(len=7), parameter :: objFolder = "../obj/"
    character(len=8), parameter :: dataFolder = "../data/"
end module rt_constants


module rt_types

    use rt_constants
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
    
end module rt_types       


module rt_global

	use rt_constants
	use rt_types
	implicit none
	real(dp), dimension(:), allocatable           :: absorbed
    real(dp), dimension(:,:), allocatable         :: vertices
    type(tetraElement), dimension(:), allocatable :: tetraData
    
end module







        
    
    
    