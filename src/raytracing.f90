! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 21.02.2014
! version: 0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use pre_process_data

    implicit none
    character(len=100)                               :: file_name
    character(len=100), dimension(1)                 :: emSurfNames
    type(rayContainer)                               :: ray
    type(tetraElement), dimension(:), allocatable    :: tetraData
    type(emissionSurface), dimension(:), allocatable :: emSurf
    real(dp), dimension(:,:), allocatable            :: vertices
    integer                                          :: npart
    
    ! user input
    file_name = "sphere.msh"
    npart = 4
    emSurfNames = ['xLow'] ! right now this requires to state the correct dimension in above declaration       
        
    ! perform pre-processing
    call read_mesh_data(file_name, emSurfNames, npart, tetraData, vertices, emSurf)
    
    ! create ray
    call CreateRay(tetraData, vertices, emSurf(1), ray)
    call TraceRay(tetraData, vertices, ray)
    
    contains 
 
    ! given an emission surface, the following subroutine perform following steps:
    ! 1. select a face on the emission surface by roulette wheel selection
    ! 2. select a point within the triangular face based on 2 random numbers 
    !    (uses triangle and uniform distribution)
    ! 3. select a random direction
    subroutine CreateRay(tetraData, vertices, ems, ray)

    ! procedure for creating starting point for ray on emission surface ems
    ! procedure chooses a random tetraeder on the surface and the creates a
    ! random point with the face on the surface
    
        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(emissionSurface), intent(in)   :: ems
        type(rayContainer), intent(out)     :: ray
        integer                             :: i
        integer, dimension(3)               :: vertIDs 
        real(dp)                            :: psi, theta, b, c, d
        real(dp), dimension(3)              :: p1, p2, p3, dir1, dir2, dir21, ds1, ndir
        real(dp), dimension(3,3)            :: M  

        ! get random tetraeder on the emission surface
        psi = myRandom(2903) ! initialize random generator
        
        ! decision whether to start at the beginning or end of the list just for speed
        if (psi > 0.5) then
            do i = size(ems%area)-1,1,-1
                if (ems%area(i) < psi) exit
            end do
            ray%tetraID = ems%elemData(i+1,1)
            ray%faceID = ems%elemData(i+1,2)
        else           
            do i = 1,size(ems%area)
                if (ems%area(i) > psi) exit
            end do
            ray%tetraID = ems%elemData(i,1)
            ray%faceID = ems%elemData(i,2)
        end if
    
        ! get vertices for the selected face on the surface
        call return_facevertIds(ray%faceID,vertIDs)  
        p1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
        p2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:)
        p3 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:)
     
        ! find longest side of triangle
        if ((norm(p2-p1) >= norm(p3-p1)) .and. (norm(p2-p1) >= norm(p3-p2))) then
            dir1 = p2-p1
            dir2 = p3-p1
            ray%point = p1     
!             write(*,*) "p1"
!             write(*,*)   
        else if ((norm(p1-p3) >= norm(p2-p1)) .and. (norm(p1-p3) >= norm(p3-p2))) then
            dir1 = p1-p3
            dir2 = p2-p3 
            ray%point = p3
!             write(*,*) "p3"
!             write(*,*)  
        else
            dir1 = p3-p2
            dir2 = p1-p2
            ray%point = p2
!             write(*,*) "p2"
!             write(*,*)
        end if
        
        ! projection of dir2 onto dir1
        dir21 = dot_product(dir1,dir2)/dot_product(dir1,dir1) * dir1
    
        ! choose random numbers based on triangle distribution
        psi = myRandom(0)
        theta = myRandom(0)

        b = norm(dir1)     ! length of dir1
        c = norm(dir21)    ! length of projection dir21
        ds1 = dir2 - dir21 ! perpendicular to d1 and going through vertex of triangle which is not on d1

        if (psi <= c/b) then
            ray%point = ray%point + sqrt(psi*b*c)*dir1/b                       ! random number along longest side (triangle distribution)
            d = sqrt(psi*b*c)*tan(acos(dot_product(dir1,dir2)/(b*norm(dir2)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d*ds1/norm(ds1)                      ! random number along perpendicular direction (uniform distribution)     
        else
            ray%point = ray%point + (b - sqrt(b*(b-c)*(1-psi)))*dir1/b                               ! random number along longest side (triangle distribution)
            d = sqrt(b*(b-c)*(1-psi)) *tan(acos(dot_product(-dir1,(dir2-dir1))/(b*norm(dir2-dir1)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d*ds1/norm(ds1)                                            ! random number along perpendicular direction (uniform distribution)        
        end if
        
        ! test whether point is indeed in triangle
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Emission Point not in chosen face!"
            stop
        end if
    
        ! choose ray direction
        ndir = cross(dir1,dir2) ! normal vector of triangle pointing outwards
        ndir = ndir/norm(ndir)
        
        ! check whether normal vector points inwards or outwards
        do i = 1,4
            if (any(i == vertIDs) .eqv. .false.) exit
        end do        
        if (dot_product(ndir, vertices(tetraData(ray%tetraID)%vertexIds(i),:) - ray%point) < 0) ndir = -ndir

        ! first rotate about about a vector in the plane
        dir1 = dir1/b
        call RotationMatrix(dir1, ndir, asin(sqrt(myRandom(0))), M)
        ray%direction = matmul(M,ndir)
        
        ! second rotation with normal vector as rotation axis
        call RotationMatrix(ndir, dir1, asin(sqrt(myRandom(0))), M)
        ray%direction = matmul(M,ray%direction)
        
        write(*,*)
        write(*,'(a12,3(1x,3e14.6))') "origin: ", ray%point
        write(*,'(a12,3(1x,3e14.6))') "direction: ", ray%direction
        write(*,*)         
        
    end subroutine CreateRay
    
    subroutine TraceRay(tetraData, vertices, ray)
    
        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        real(dp)                            :: kappa, sigma ! absorption and scattering coefficients (must be provided from somewhere)
        real(dp)                            :: lAbs, lScat, alpha, temp
        integer                             :: i, newFace
        integer, dimension(3)               :: vertIDs 
        real(dp), dimension(3)              :: rp, v1, v2, nsf
        logical                             :: test
        
        ! calculate length of initial ray
        kappa = 1.0_dp ! simple assumption so far
        sigma = 1.0_dp ! simple assumption so far
        
        lAbs = 1.0_dp/kappa*log(1/myRandom(0))
        lScat = 1.0_dp/sigma*log(1/myRandom(0))
        
        if (lAbs < lScat) then
            
            ! perform absorption
        else
            
            ! perform scattering
        
        end if
        
        write(*,*) lAbs
        write(*,*) lScat
        
        ! below: current work ....
        
        write(*,*)
        write(*,*) ray%faceID
        write(*,*) ray%tetraID
        
        ! trace ray
        alpha = 100.0_dp  ! some large initial value for alpha
        do i = 1,4
        
            if (i == ray%faceID) cycle  ! is face where ray is emitted
            
            ! 1 vertex and 2 vectors of current face
            call return_facevertIds(i,vertIDs)
            rp = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
            v1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:) - rp  
            v2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:) - rp
            
            ! normal vector of current face  
            nsf = cross(v1,v2)
            nsf = nsf/norm(nsf) 

            ! calculate length until intersection
            temp = (dot_product(nsf,rp) - dot_product(nsf,ray%point))/dot_product(nsf, ray%direction)
            if (temp < 0) cycle
            ! check if it is the shortest length
            if (temp < alpha) then
                alpha = temp
                newFace = i
            end if
            
            write(*,*) alpha
            
        end do
        
        ray%point = ray%point + alpha*ray%direction
        
        ! get tetraeder which shares the same face
        temp = tetraData(ray%tetraID)%neighbors(newFace,1)
        ray%faceID = tetraData(ray%tetraID)%neighbors(newFace,2)
        ray%tetraID = temp    
        
        write(*,*)
        write(*,*) ray%faceID
        write(*,*) ray%tetraID
        
        call return_facevertIds(ray%faceID, vertIDs)
        rp = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
        v1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:)   
        v2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:)
        test = PointInside(rp,v1,v2,ray%point)
        write(*,*) test
        
    end subroutine TraceRay
    
    ! test whether a point is inside a triangle or not
    function PointInside(p1,p2,p3,tp)
        
        real(dp), dimension(3), intent(in) :: p1,p2,p3,tp
        logical                            :: PointInside
        real(dp), dimension(3)             :: na,nx,ny,nz,bc
        real(dp)                           :: dpna
        
        na = cross(p2-p1,p3-p1)
        nx = cross(p3-p2,tp-p2)
        ny = cross(p1-p3,tp-p3)
        nz = cross(p2-p1,tp-p1)
        dpna = dot_product(na,na)
        bc(1) = dot_product(na,nx)/dpna
        bc(2) = dot_product(na,ny)/dpna
        bc(3) = dot_product(na,nz)/dpna
        PointInside = (all(bc > 0e0_dp) .and. (abs(sum(bc) - 1) <= 1e-14_dp))
        
!         write(*,*) bc(1)
!         write(*,*) bc(2)
!         write(*,*) bc(3)
!         write(*,*) 1-sum(bc)
    
    end function
  
  
end program raytracing 