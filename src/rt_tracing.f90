! Raytracing tracing module
! author: Steffen Finck
! contact: steffen.finck@fhv.at

module tracing

    use rt_funcs
    use math_funs
    use helper_functions
    
    implicit none
    
    contains
    
    subroutine start_tracing(tetraData, vertices, ems, file_name, nrays)
    
		type(tetraElement), intent(inout)      :: tetraData(:)
        real(dp), intent(inout)                :: vertices(:,:)
        type(emissionSurface), intent(inout)   :: ems(:)
        character(len=*), intent(in)           :: file_name 
        integer, intent(in)                    :: nrays
        type(rayContainer)                     :: ray
        real(dp) :: t1, t2, dummy
        integer :: io_error, write_error, k
        character(len=100) :: resFname, leaveFname
        
        ! initialize random generator
        dummy = myRandom(2903)
        
        ! create file for output
        ! file containing absorbing location
        resFname = objFolder//trim(file_name)//".res"
        open(unit=81, file=resFname, action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating file for results 1",81)
        write(81,'(a8,4(1x,a14))',iostat=write_error) "tetraID", "x", "y", "z", "absorbed"
        call check_io_error(write_error,"writing header results 1",81)
        close(unit=81, iostat=io_error)
        call check_io_error(io_error,"closing file for results 1",81)
        
        ! file containing leaving enclosure location
        leaveFname = objFolder//trim(file_name)//"-leaving.res"
        open(unit=83, file=leaveFname, action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating file for results 2",83)
        write(83,'(7(1x,a14))',iostat=write_error) "loc x", "loc y", "loc z", "dir x", "dir y", "dir z", "power" 
        call check_io_error(write_error,"writing header results 2",83)
        close(unit=83, iostat=io_error)
        call check_io_error(io_error,"closing file for results 2",83)
        
        ! additional debugging files
        open(unit=82, file=objFolder//trim(file_name)//"-initrayloc.res", action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating debug-file 1",82)
        write(82,'(6(1x,a14))',iostat=write_error) "x", "y", "z", "dx", "dy","dz"
        call check_io_error(write_error,"writing debug-file 1 header",82)
        close(unit=82, iostat=io_error)
        call check_io_error(io_error,"closing debug-file 1",82)
        
        ! do raytracing            
        ! in case of several emission surface a way of selecting a surface must be chosen     
	    call cpu_time(t1)
	    do k = 1,nrays
	        call CreateRay(tetraData, vertices, ems(1), ray, file_name)
	        call TraceRay(tetraData, vertices, ray, resFname, leaveFname)
	    end do
	    call cpu_time(t2)
    
	    write(*,*) "run time raytracing: ", t2-t1
    
    end subroutine start_tracing
    
    ! given an emission surface, the following subroutine performs the following steps:
    ! 1. select a face on the emission surface by roulette wheel selection
    ! 2. select a point within the triangular face based on 2 random numbers 
    !    (uses triangle and uniform distribution)
    ! 3. select a random direction    
    subroutine CreateRay(tetraData, vertices, ems, ray, file_name)

        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(emissionSurface), intent(in)   :: ems
        type(rayContainer), intent(out)     :: ray
        character(len=*), intent(in)        :: file_name
        integer  :: i
        integer, dimension(1) :: id
        integer, dimension(3) :: vertIDs 
        real(dp) :: psi, theta, b, c, d
        real(dp), dimension(3) :: p1, p2, p3, dir1, dir2, dir21, ds1, ndir
        real(dp), dimension(3,3) :: M1, M2  

        ! get random tetraeder on the emission surface
        psi = myRandom(0) 
        
        ! decision whether to start at the beginning or end of the list (just for speed)
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
        call return_coords(tetraData(ray%tetraID), vertices, vertIDs, p1, p2, p3)
     
        ! find longest side of triangle
        id = maxloc([norm(p2-p1),norm(p3-p1),norm(p2-p3)])        
        select case (id(1))
	        case(1)
		        dir1 = p2-p1
	            dir2 = p3-p1
	            ray%point = p1
	        case(2)
		        dir1 = p1-p3
                dir2 = p2-p3 
                ray%point = p3
            case default
	            dir1 = p3-p2
                dir2 = p1-p2
                ray%point = p2
        end select
        
        ! projection of dir2 onto dir1
        dir21 = dot_product(dir1,dir2)/dot_product(dir1,dir1) * dir1
    
        ! choose random numbers based on triangle distribution
        psi = myRandom(0)
        theta = myRandom(0)
 
        ! parameters of triangle distribution
        b = norm(dir1)     ! length of dir1
        c = norm(dir21)    ! length of projection dir21
        ds1 = dir2 - dir21 ! perpendicular to d1 and going through vertex of triangle which is not on d1
                           ! dir21 + ds1 = dir2 
        
        ! choose length depending on cdf of triangle distribution
        if (psi <= c/b) then
            ray%point = ray%point + sqrt(psi*b*c)*dir1/b                       ! random number along longest side (triangle distribution)
            d = sqrt(psi*b*c)*tan(acos(dot_product(dir1,dir2)/(b*norm(dir2)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d/norm(ds1)*ds1                      ! random number along perpendicular direction (uniform distribution)     
        else
            ray%point = ray%point + (b - sqrt(b*(b-c)*(1-psi)))*dir1/b                               ! random number along longest side (triangle distribution)
            d = sqrt(b*(b-c)*(1-psi)) *tan(acos(dot_product(-dir1,(dir2-dir1))/(b*norm(dir2-dir1)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d/norm(ds1)*ds1                                            ! random number along perpendicular direction (uniform distribution)        
        end if
        
        ! test whether point is indeed in triangle
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Emission Point not in chosen face!"
            stop
        end if
    
        ! choose ray direction
        ndir = cross(dir1,dir2) ! normal vector of triangle face
        ndir = ndir/norm(ndir)
        
        ! check whether normal vector points inwards or outwards
        do i = 1,4
            if (any(i == vertIDs) .eqv. .false.) exit
        end do        
        if (dot_product(ndir, vertices(tetraData(ray%tetraID)%vertexIds(i),:) - ray%point) < 0) ndir = -ndir

        ! rotate normal vector into direction defined by psi and theta
        dir1 = dir1/b
        theta = asin(sqrt(myRandom(0)))
        psi = 2.0_dp*pi*myRandom(0) 	         
        ray%direction = sin(theta)*(cos(psi)*dir1 + sin(psi)*cross(ndir, dir1)) + cos(theta)*ndir
        
!         dir1 = dir1/b
!         call RotationMatrix(dir1, ndir, , M1)
!     
!         ! second rotation with normal vector as rotation axis
!         call RotationMatrix(ndir, dir1, , M2)
!         
!         ! get final direction
!         ray%direction = matmul(M2,matmul(M1,ndir))
        
        ! just for checking (could be commented)
        open(unit=83, file=objFolder//trim(file_name)//"-initrayloc.res", action='write', position='append')  
        write(83,'(6(1x,e14.6))') ray%point, ray%direction
        close(unit=83)         
        
    end subroutine CreateRay
    
    ! main raytracing routine
    subroutine TraceRay(tetraData, vertices, ray, resFname, leaveFname)
    
        type(tetraElement), intent(inout)   :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        character(len=*), intent(in)        :: resFname, leaveFname
        real(dp) :: kappa, sigma, beta, omega, length
        integer :: flag, cface
        real(dp), dimension(3) :: ipoint 
        type(tetraElement) :: tetra
        
        ! properties of the medium
        kappa = 10.0_dp
        sigma = 0.75_dp
        
        ! calculate length of initial ray
        beta = kappa + sigma        
        omega = sigma/beta
        
        ! trace ray until the energy is below a treshold  
        do while (ray%power > 0.00001_dp)
        
	        length = 1.0_dp/beta*log(1/myRandom(0))
	        ipoint = ray%point + length*ray%direction
	        flag = 1
            
            ! trace path 
            do while (ray%length < length)                
                ! current tetra
                tetra = tetraData(ray%tetraID)       
		        cface = ray%faceID ! to avoid association behavior in function call below
		        
		        if (cface <0 ) stop
		        
		        ! find next face within same tetra
		        call FindNextFace(vertices,tetra,ray,cface)
		         
		        ! check if point is on a boundary surface
	            if (tetra%neighbors(ray%faceID,2) < 0) then
	                call BoundaryHandling(ray, tetra, vertices, leaveFname)
                    flag = 0
                    exit
			    end if
			    
			    ! check if mean-free path length is reached
			    if (ray%length >= length) exit
			    
			    ! update remaining raycontainer values
			    ray%tetraID = tetra%neighbors(ray%faceID,1) ! neighbouring tetra
			    ray%faceID = tetra%neighbors(ray%faceID,2)  ! face in neighbouring tetra
            end do
            
            ! ray properties are updated and no scattering or absorption happens
            if (flag < 1) cycle
            
            ! set ray point to point where interaction happens
            ray%point = ipoint
            
            !absorption
            call RayAbsorbing(tetra,ray%tetraID,ray%point,(1.0_dp-omega)*ray%power,resFname)
	        ray%power = ray%power*omega
	        
	        ! scattering
	        call RayScatter(ray,tetra,vertices)
	        ! if scattering to a boundary happens
	        if (tetra%neighbors(ray%faceID,2) < 0) then
		        call BoundaryHandling(ray, tetra, vertices, leaveFname)
		        cycle
		    end if
		    
	        ! update ray container
	        ray%tetraID = tetra%neighbors(ray%faceID,1) ! neighbouring tetra
	        ray%faceID = tetra%neighbors(ray%faceID,2)  ! face in neighbouring tetra
	    
	        ! update ray properties    
		    ray%length = 0.0_dp
            
        end do 
        
        ! in case a small amount of power remains
        if (ray%power > 0) call RayAbsorbing(tetraData(ray%tetraID),ray%tetraID,ray%point,ray%power,resFname)
        
    end subroutine TraceRay
    
    
    ! find next next face along path in same tetra
    subroutine FindNextFace(vertices,tetra,ray,faceOld)
        
        type(tetraElement), intent(in)      :: tetra
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        integer, intent(in)                 :: faceOld
        real(dp), dimension(3) :: p1, p2, p3, nsf
        integer, dimension(3)  :: vertIDs
        real(dp) :: alpha, minalpha
	    integer :: f
                        
        ! get new point on face of tetraeder
        minalpha = 1000.0_dp
        do f = 1,4
            
            ! if emitted from, cycle if it is the same face
            if (f == faceOld) cycle  
            
            ! 1 vertex and 2 vectors of current face
            call return_facevertIds(f,vertIDs)
            call return_coords(tetra, vertices, vertIDs, p1, p2, p3)
            p2 = p2 - p1 ! make a vector with origin at p1  
            p3 = p3 - p1 ! make a vector with origin at p1
            
            ! normal vector of current face  
            nsf = cross(p2,p3)
            nsf = nsf/norm(nsf) 

            ! calculate length until intersection
            alpha = (dot_product(nsf,p1) - dot_product(nsf,ray%point))/dot_product(nsf, ray%direction)
            if (alpha < 0) cycle ! length can not be negative
            
            ! check if temp is the shortest length so far
            ! if true, current face is the face the ray will intersect
            if (alpha < minalpha) then
                minalpha = alpha
                ray%faceID = f
            end if
            
        end do
        
        ! update some raycontainer elements
        ray%point = ray%point + minalpha*ray%direction
        ray%length = ray%length + minalpha            
    
        ! test whether point is indeed in triangle
        call return_facevertIds(ray%faceID,vertIDs)  
        call return_coords(tetra, vertices, vertIDs, p1, p2, p3)
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Waypoint not in face!"
            write(*,*) "rp:", ray%point
            write(*,*) "rt:", ray%tetraID
            write(*,*) "rf:", ray%faceID
            write(*,*) "rd:", ray%direction
            write(*,*) "alpha:", minalpha
            write(*,*) "p1:", p1
            write(*,*) "p2:", p2
            write(*,*) "p3:", p3
            write(*,*) "fold:", faceOld
            stop
        end if
    
    end subroutine FindNextFace
        
        
    ! write out location and intensity of ray absorbed
    subroutine RayAbsorbing(tetra, id, point, power, resFname)
    
	    type(tetraElement), intent(inout)  :: tetra
	    integer, intent(in)                :: id
	    real(dp), dimension(3), intent(in) :: point
	    real(dp), intent(in)               :: power
	    character(len=*), intent(in)       :: resFname
	    
	    ! do some shape function magic
	    
	    open(unit=84, file=resFname, action='write', position='append')  
	    write(84,'(1x,i8,1x,3(e14.6,1x),e14.6)') id, point, power
        close(unit=84)
        
    end subroutine RayAbsorbing
    
    
    ! perform scattering (so far isotropic only)
    subroutine RayScatter(ray, tetra, vertices)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    real(dp), intent(in)              :: vertices(:,:)
	    real(dp) :: psi, theta
	    real(dp), dimension(3) :: v1
	    
	    ! isotropic case
	    theta = acos(1.0_dp-2.0_dp*myRandom(0))
	    psi = 2.0_dp*pi*myRandom(0)
		
! 		write(*,*) "theta:", theta
! 		write(*,*) "psi:", psi
		
		! determine new direction
		! get 2 perpendicular vectors, such that the old direction is the normal vector
		! for the plane spanned by the 2 new vectors
		v1 = cross(eoshift(ray%direction,1),ray%direction)
		v1 = v1/norm(v1)	         
        ray%direction = sin(theta)*(cos(psi)*v1 + sin(psi)*cross(ray%direction, v1)) + cos(theta)*ray%direction
        
        ! find intersection point with face of current tetraeder
        call FindNextFace(vertices,tetra,ray,-1)
        
        ! some test of sanity?
        
    end subroutine RayScatter
        
        
    ! calculate refraction with Snell's law
    subroutine SnellsLaw(incident, n1, n2, nsf, reflect, refract)
        
        real(dp), intent(in)                  :: n1, n2   ! refractive indices
        real(dp), dimension(3), intent(in)    :: incident ! incident direction 
        real(dp), dimension(3), intent(inout) :: nsf      ! normal vector
        real(dp), dimension(3), intent(out)   :: refract, reflect ! refraction and reflection direction
        real(dp) :: cosAngle, ratio
        
        ! cosine of angle between surface normal and incident direction
        ! must be positive (else normal is pointing outwards and must be used as negative)
        cosAngle = -dot_product(incident,nsf)
        if (cosAngle < 0.0_dp) then
            cosAngle = -cosAngle
            nsf = -nsf
        end if
        
        ! direction of reflected ray
        reflect = incident + 2*cosAngle*nsf
        
        ! ratio of refractive indices
        ratio = n1/n2
        
        ! check whether critical angle is exceeded 
        if (ratio*sin(acos(cosAngle)) > 1) then
	        refract = [0.0_dp,0.0_dp,0.0_dp]
	    else
	        refract = ratio*incident + (ratio*cosAngle - sqrt(1.0_dp - ratio**2*(1.0_dp-cosAngle**2)))*nsf
        end if
        
    end subroutine SnellsLaw
    
    ! handle boundary
    ! parameter trans specifies transmissivity and reflectivity (1-trans)
    subroutine BoundaryHandling(ray, tetra, vertices, leaveFname)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    real(dp), intent(in)              :: vertices(:,:)
	    character(len=*), intent(in)      :: leaveFname
	    real(dp), dimension(3) :: nsf, point, reflect, refract 
	    
	    ! some sanity check
	    if (ray%faceID < 0 ) then
		    write(*,*) " negative face id, where are you coming from?"
		    stop
		end if
	    
	    ! get surface normal of current face
	    call return_surfNormal(tetra, ray%faceID, vertices, nsf, point)
	    ! reflection and refraction  
	    call SnellsLaw(ray%direction, 1.9_dp, 1.5_dp, nsf, reflect, refract)
! 	    call SnellsLaw(refract, 1.5_dp, 1.0_dp, nsf, reflect, refract)
	    
	    ray%direction = reflect
	    if (count(refract == 0.0_dp) < 3) then		    
		    
		    ! write out data 
		    open(unit=85, file=leaveFname, action='write', position='append')  
		    write(85,'(7(1x,e14.6))') ray%point, refract, ray%power*trans
	        close(unit=85)
		    
		    ! update ray properties
	        ray%power = ray%power*(1-trans)
	        
	    end if
	    
	    ! since current face id is negative, reset it
!         ray%faceID = get_facenumber(ray%faceID, tetra)
        ray%length = 0.0_dp
        
!         write(*,*) "rp", ray%point
!         write(*,*) "rd", ray%direction
	           
	end subroutine BoundaryHandling
	
end module tracing