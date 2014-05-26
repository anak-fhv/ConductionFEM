! Raytracing tracing module
! author: Steffen Finck
! contact: steffen.finck@fhv.at

module tracing

	use rt_elements
    use rt_helper
    
    implicit none
    
    contains
    
    ! wrapper routine and does perform initialization of fields, output files, etc.
    subroutine start_tracing
    
        type(rayContainer) :: ray
        real(dp)           :: t1, t2
        integer            :: io_error, write_error, k, alloc_status
        character(len=100) :: leaveFname, d_file1, d_file2, tmp 
        
        ! create file for output
        ! file containing leaving enclosure location
        leaveFname = objFolder//trim(data_fname)//"-leaving.res"
        open(unit=83, file=leaveFname, action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating file for results 2",83)
        write(83,'(8(1x,a14))',iostat=write_error) "loc x", "loc y", "loc z", "dir x", "dir y", "dir z", "power", "wavelength" 
        call check_io_error(write_error,"writing header results 2",83)
        close(unit=83, iostat=io_error)
        call check_io_error(io_error,"closing file for results 2",83)
        
        ! additional debugging files
        ! distribution of initial locations
        d_file1 = objFolder//trim(data_fname)//"-initrayloc.res" 
        open(unit=82, file=d_file1, action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating debug-file 1",82)
        write(82,'(7(1x,a14))',iostat=write_error) "x", "y", "z", "dx", "dy","dz", "power"
        call check_io_error(write_error,"writing debug-file 1 header",82)
        close(unit=82, iostat=io_error)
        call check_io_error(io_error,"closing debug-file 1",82)
             
        ! initialize vector for absorption
        allocate(absorbed(size(vertices,1)), stat=alloc_status)
        call check_alloc_error(alloc_status, "absorbed vector")
        absorbed = 0.0_dp
        
        ! do raytracing          
        call cpu_time(t1)
	    do k = 1,nrays
		    
		    if (k <= nRayPaths) then
			    ! file for tracing single ray (should be commented if large number of rays are considered)
			    write(tmp,'(i5.5)') k
			    d_file2 = objFolder//trim(data_fname)//"-raydata"//trim(tmp)//".res" 
		        open(unit=81, file=d_file2, action='write', status='replace', iostat=io_error)       
		        call check_io_error(io_error,"creating ray-tracing file",81)
		        write(81,'(7(1x,a14))',iostat=write_error) "x", "y", "z", "dx", "dy","dz","power" 
		        call check_io_error(write_error,"writing ray-tracing file header",81)
	            close(unit=81, iostat=io_error)
		        call check_io_error(io_error,"closing ray-tracing file",81)
		    end if
		    
		    ! create ray
	        call CreateRay(emSurf(emsIDfun()), ray, d_file1)

	        ! assign energy to a ray and trace it
            if (k <= nRayPaths) call WriteRayData(ray, d_file2)
	        call TraceRay(ray, leaveFname, d_file2, k<=nRayPaths)
	        
	    end do
	    call cpu_time(t2)
    
	    write(*,*) "run time raytracing: ", t2-t1
    
    end subroutine start_tracing
    
    
    ! given an emission surface, the following subroutine performs the following steps:
    ! 1. select a face on the emission surface by roulette wheel selection
    ! 2. select a point within the triangular face based on 2 random numbers 
    !    (uses triangle and uniform distribution)
    ! 3. select a random direction    
    subroutine CreateRay(ems, ray, fname)

        type(emissionSurface), intent(in) :: ems
        type(rayContainer), intent(out)   :: ray
        character(len=*), intent(in)      :: fname
        integer                           :: i
        integer, dimension(1)             :: id
        integer, dimension(3)             :: vertIDs 
        real(dp)                          :: psi, theta, area, tc1, tc2, cTemperature
        real(dp), dimension(3)            :: p1, p2, p3, dir1, ndir 
        type(tetraElement)                :: tetra
        
        ! decision whether to start at the beginning or end of the list (just for speed)
        psi = myRandom(0)
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
	    tetra = tetraData(ray%tetraID)
    
        ! get vertices for the selected face on the surface
        call return_facevertIds(ray%faceID,vertIDs)  
        call return_coords(tetra, vertIDs, p1, p2, p3)
        
        ! get random value for two tetra-coordinates
        tc1 = 1.0_dp - sqrt(1-myRandom(0))
        tc2 = (1.0_dp-tc1)*myRandom(0)
        
        ! get points from triangle coordinates
        ray%point = tc1*p1 + tc2*p2 + (1.0_dp - tc1 - tc2)*p3 
        
        ! test whether point is indeed in triangle
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Emission Point not in chosen face!"
            stop
        end if
    
        ! choose ray direction
!         ndir = cross(dir1,dir2)  ! normal vector of triangle face
        ndir = cross(p2-p1, p3-p1)
        area = 0.5_dp*norm(ndir) ! area of triangle
        ndir = ndir/norm(ndir)
        
        ! check whether normal vector points inwards or outwards
        if (dot_product(ndir, vertices(tetra%vertexIds(10-sum(vertIDs)),:) - ray%point) < 0) ndir = -ndir

        ! rotate normal vector into direction defined by psi and theta
!         dir1 = dir1/b
        dir1 = (p2-p1)/norm(p2-p1)
        theta = asin(sqrt(myRandom(0)))
        psi = 2.0_dp*pi*myRandom(0) 	         
        ray%direction = sin(theta)*(cos(psi)*dir1 + sin(psi)*cross(ndir, dir1)) + cos(theta)*ndir
        
        ! power of the ray
        ! ray%power = RayPowerFun(tc1*tData(tetra%vertexIds(vertIDs(1))) + tc2*tData(tetra%vertexIds(vertIDs(2))) + (1-tc1-tc2)*tData(tetra%vertexIds(vertIDs(3))),area)
        cTemperature = tc1*temperature(tetra%vertexIds(vertIDs(1))) + tc2*temperature(tetra%vertexIds(vertIDs(2))) + (1.0_dp - tc1- tc2)*temperature(tetra%vertexIds(vertIDs(3)))
        !write(*,*) cTemperature
        if (count(ems%name == ignoredSurfaces) == 0) then
	        ray%power = RayPowerFun(cTemperature,area, 1.0_dp)
	    else
		    ray%power = RayPowerFun(cTemperature,area, alpha)
		end if
        
        Etotal = Etotal  +  ray%power
        
        ! just for checking (could be commented)
        open(unit=83, file=fname, action='write', position='append')  
        write(83,'(7(1x,e14.6))') ray%point, ray%direction, ray%power
        close(unit=83)         
        
    end subroutine CreateRay
    
    
    ! main raytracing routine
    subroutine TraceRay(ray, leaveFname, rtfname, writeflag)
    
        type(rayContainer), intent(inout) :: ray
        character(len=*), intent(in)      :: leaveFname, rtfname
        logical, intent(in)               :: writeflag
        real(dp)                          :: beta, omega, length, tmp
        integer                           :: cface 
        type(tetraElement)                :: tetra
        
        ! calculate length of initial ray
        ! kappa and sigma are defined by module rt_properties
        ! only necessary in case of participating media
        beta = kappa + sigma        
        omega = sigma/beta
        
        ! trace ray until the energy is below a treshold  
        tmp = ray%power  ! initial ray power
        length = GetPathLength() ! path length (only meaningful for participating media)
        do while (ray%power > 0.0001_dp*tmp)
        
	        ! current tetra
            tetra = tetraData(ray%tetraID)       
		    cface = ray%faceID 
		        		        
		    ! find next face within same tetra
		    call FindNextFace(tetra,ray,cface)
            if (writeflag) call WriteRayData(ray, rtfname)
                
            ! check if point is on a boundary surface
	        if (tetra%neighbors(ray%faceID,2) < 0) then
	        
	            call BoundaryHandling(ray, tetra, leaveFname)
	            
	            ! write ray path data	            
	            if (writeflag) call WriteRayData(ray, rtfname)
	            
	            ! reset ray%length and get new mean free path
                ray%length = 0.0_dp
				length = GetPathLength()
				
            elseif ((RT_setup .eq. 'led') .and. (ray%length >= length) ) then 
			    
			    ! setup specific for 'led'    
				! check if mean-free path length is reached
				
			    ! set ray point to point where interaction happens
			    ray%point = ray%point + (length - ray%length)*ray%direction
            
			    !absorption
			    call RayAbsorbing(ray, tetra, 1.0_dp-omega)
	        
				! scattering
			    call RayScatter(ray,tetra)
				! if scattering to a boundary happens
				if (tetra%neighbors(ray%faceID,2) < 0) then
				    call BoundaryHandling(ray, tetra, leaveFname)
				end if
						
				! write point on path
			    if (writeflag) call WriteRayData(ray, rtfname)
					    
				! reset ray%length and get new mean free path length
				ray%length = 0.0_dp
				length = GetPathLength()	    
                
			elseif ((RT_setup .eq. 'tomo') .and. (tetra%domain /= tetraData(tetra%neighbors(ray%faceID,1))%domain)) then  
			
				! setup specific for 'tomo'  
				! check if domain will change
				
			    call DomainChange(ray, tetra)
	            if (writeflag) call WriteRayData(ray, rtfname)
				
			else     ! nothing happens just go on with tracing
			
				! update ray identifier values
			    ray%tetraID = tetra%neighbors(ray%faceID,1) ! neighbouring tetra
			    ray%faceID = tetra%neighbors(ray%faceID,2)  ! face in neighbouring tetra
				    
		    end if 			    
			    			    			    
        end do                        
		    
        ! in case a small amount of power remains, absorb ray completely
        if (ray%power > 0.0_dp) then
	        call RayAbsorbing(ray, tetra, 1.0_dp)
	        if (writeflag) call WriteRayData(ray, rtfname)    
	    end if        
        
    end subroutine TraceRay
    
    
    ! find next next face along path in same tetra
    subroutine FindNextFace(tetra,ray,faceOld)
        
        type(tetraElement), intent(in)    :: tetra
        type(rayContainer), intent(inout) :: ray
        integer, intent(in)               :: faceOld
        real(dp), dimension(3)            :: p1, p2, p3, nsf
        integer, dimension(3)             :: vertIDs
        real(dp)                          :: length, minlength
	    integer                           :: f
                        
        ! get new point on face of tetraeder
        minlength = 1000.0_dp
        do f = 1,4
            
            ! if emitted from, cycle if it is the same face
            if (f == faceOld) cycle  
            
            call return_surfNormal(tetra, f, nsf, p1)
            
!             ! 1 vertex and 2 vectors of current face
!             call return_facevertIds(f,vertIDs)
!             call return_coords(tetra, vertIDs, p1, p2, p3)
!             p2 = p2 - p1 ! make a vector with origin at p1  
!             p3 = p3 - p1 ! make a vector with origin at p1
!             
!             ! normal vector of current face  
!             nsf = cross(p2,p3)
!             nsf = nsf/norm(nsf) 

            ! calculate length until intersection
            length = (dot_product(nsf,p1) - dot_product(nsf,ray%point))/dot_product(nsf, ray%direction)
            if (length < 0) cycle ! length can not be negative
            
            ! check if temp is the shortest length so far
            ! if true, current face is the face the ray will intersect
            if (length < minlength) then
                minlength = length
                ray%faceID = f
            end if
            
        end do
        
        ! update some raycontainer elements
        ray%point = ray%point + minlength*ray%direction
        ray%length = ray%length + minlength            
    
        ! test whether point is indeed in triangle
        call return_facevertIds(ray%faceID,vertIDs)  
        call return_coords(tetra, vertIDs, p1, p2, p3)
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Waypoint not in face!"
            write(*,*) "rp:", ray%point
            write(*,*) "rt:", ray%tetraID
            write(*,*) "rf:", ray%faceID
            write(*,*) "rd:", ray%direction
            write(*,*) "alpha:", minlength
            write(*,*) "p1:", p1
            write(*,*) "p2:", p2
            write(*,*) "p3:", p3
            write(*,*) "fold:", faceOld
            
            call return_facevertIds(faceOld,vertIDs)  
            call return_coords(tetra, vertIDs, p1, p2, p3)
            
            write(*,*) "p4:", p1
            write(*,*) "p5:", p2
            write(*,*) "p6:", p3
            
            stop
        end if
    
    end subroutine FindNextFace
        
        
    ! write out location and intensity of ray absorbed
    subroutine RayAbsorbing(ray, tetra, frac)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    real(dp), intent(in)              :: frac
	    real(dp)                          :: t1, t2, t3
	    
	    ! do some shape function magic
	    call Cartesian2Tetra(ray%point, tetra, t1, t2, t3)
	    
	    ! save absorption values
	    absorbed(tetra%vertexIds(1)) = absorbed(tetra%vertexIds(1)) + t1*frac*ray%power
        absorbed(tetra%vertexIds(2)) = absorbed(tetra%vertexIds(2)) + t2*frac*ray%power
        absorbed(tetra%vertexIds(3)) = absorbed(tetra%vertexIds(3)) + t3*frac*ray%power
        absorbed(tetra%vertexIds(4)) = absorbed(tetra%vertexIds(4)) + max((1-t1-t2-t3),0.0_dp)*frac*ray%power
        
        ! update raypower
        ray%power = (1-frac)*ray%power
!         write(*,'(2(e14.6,1x),4(i8,1x))') ray%power, frac, tetra%vertexIds(1), tetra%vertexIds(2),tetra%vertexIds(3),tetra%vertexIds(4)
        
    end subroutine RayAbsorbing
    
    
    ! perform scattering (so far isotropic only)
    subroutine RayScatter(ray, tetra)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    real(dp)                          :: psi, theta
	    real(dp), dimension(3)            :: v1
	    
	    ! get random value for angle from phase function
	    ! phase function is defined in module rt_properties
	    call PhaseFunction(theta, psi)
		
! 		write(*,*) "theta:", theta
! 		write(*,*) "psi:", psi
		
		! determine new direction
		! get 2 perpendicular vectors, such that the old direction is the normal vector
		! for the plane spanned by the 2 new vectors
		v1 = cross(eoshift(ray%direction,1),ray%direction) !TODO: this seems like a "good" hack
		v1 = v1/norm(v1)	         
        ray%direction = sin(theta)*(cos(psi)*v1 + sin(psi)*cross(ray%direction, v1)) + cos(theta)*ray%direction
        
        ! find intersection point with face of current tetraeder
        call FindNextFace(tetra,ray,-1)
        
        ! output shoul be point on a tetraeder face with a direction inwards the tetra
        ray%tetraID = tetra%neighbors(ray%faceID,1) ! neighbouring tetra
	    ray%faceID = tetra%neighbors(ray%faceID,2)  ! face in neighbouring tetra
        
        ! some test of sanity?
        
    end subroutine RayScatter

    
    ! handle boundary
    ! either totally reflects or transmits ray depending on incoming direction
    ! and comparision with random number
    subroutine BoundaryHandling(ray, tetra, leaveFname)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    character(len=*), intent(in)      :: leaveFname
	    real(dp), dimension(3)            :: nsf
	    real(dp)                          :: cosAngle, rho, theta1, theta2, ratio
	    
	    ! TODO: 1. check on which boundary the point leaves the enclosure,
	    !         this might requires different refraction number
	    !       2. perform multiple refraction? (could be done by post-processing)
	    
	    ! some sanity check
	    if (ray%faceID < 0 ) then
		    write(*,*) " negative face id, where are you coming from?"
		    stop
		end if
		
		if (RT_setup .eq. 'tomo') then
! 		    write(*,*) 
			if (count(-tetra%neighbors(ray%faceID,2) == emSurf%originalID) == 0) then
				! one of the bounding box boundaries
				! here specular reflection occurs
				call IncidentAngle(tetra, ray, cosAngle, nsf)
				ray%direction = ray%direction + 2*cosAngle*nsf
		    else
			    ! one of the boundary emission surfaces
			    ! here complete absortption occurs
				call RayAbsorbing(ray, tetra, 1.0_dp)
		    end if
			return
		end if
		
		! refraction indices (should come from outside)
		! refraction indices are defined by module rt_properties
		if (tetra%neighbors(ray%faceID,2) < 0)  then
			ratio = refracIndices(tetra%domain+1)/refracIndices(1)
	    else
		    ratio = refracIndices(tetra%domain+1)/refracIndices(tetraData(tetra%neighbors(ray%faceID,1))%domain + 1)
		end if
	    
	    ! get incident angle
	    call IncidentAngle(tetra, ray, cosAngle, nsf)
        
        ! check for critical angle  
        if (ratio*sin(acos(cosAngle)) > 1) then
            ! total reflection
	        ray%direction = ray%direction + 2*cosAngle*nsf 
	    else
	        ! angles from Snell's law
		    theta1 = acos(cosAngle)
		    theta2 = asin(ratio*sin(theta1))
	        
	        ! Fresnel's relation
		    rho = 0.5_dp*(tan(theta1-theta2)**2/tan(theta1+theta2)**2 + sin(theta1-theta2)**2/sin(theta1+theta2)**2)
	        
	        ! decide wheter everything is transmitted or reflected
	        if (myRandom(0) < rho) then
		        ray%direction = ray%direction + 2*cosAngle*nsf  ! ray is reflected
		    else
		        ! ray is transmitted
			    ray%direction = ratio*ray%direction + (ratio*cosAngle - sqrt(1-ratio**2*(1-cosAngle**2)))*nsf 
	            
	            ! write out ray data
	            open(unit=85, file=leaveFname, action='write', position='append')  
		        write(85,'(8(1x,e14.6))') ray%point, ray%direction, ray%power, ray%wavelength
	            close(unit=85)
	            
	            ! update ray container
	            ray%power = 0.0_dp
	        end if
	        
        end if
	           
	end subroutine BoundaryHandling
	
	
	! actions and changes in case a change of domain occurs
	subroutine DomainChange(ray, tetra)

		type(rayContainer), intent(inout) :: ray
		type(tetraElement), intent(in)    :: tetra
		real(dp), dimension(3)            :: nsf, p1, dir1
		real(dp)                          :: cosAngle, theta, psi
		integer, dimension(3)             :: vertIDs
		
! 		! specular case	
! 		! get incident angle
! 		call IncidentAngle(tetra, ray, cosAngle, nsf)
!         
!         ! compare with random if reflection or absorption occurs
!         if (myRandom(0) <= 1.0_dp - 1.5_dp*alpha*cosAngle) then
! 	        ! ray is reflected
! 	        ray%direction = ray%direction + 2*cosAngle*nsf  
! 	    else
! 		    ! ray is absorbed
! 		    ! TODO: case of partial absorption
! 		    call RayAbsorbing(ray, tetra, 1.0_dp)
! 		end if
		
		! diffuse case
	    if (myRandom(0) <= 1.0_dp - alpha) then
		    ! diffuse reflection
		    call return_surfNormal(tetra, ray%faceID, nsf, p1)
		    dir1 = (ray%point-p1)/norm(ray%point-p1) ! a vector in the face
            theta = asin(sqrt(myRandom(0)))
            psi = 2.0_dp*pi*myRandom(0) 	         
            ray%direction = sin(theta)*(cos(psi)*dir1 + sin(psi)*cross(nsf, dir1)) + cos(theta)*nsf
		else
		    ! absorption
		    ! TODO: case of partial absorption
		    call RayAbsorbing(ray, tetra, 1.0_dp)
		end if
		
	end subroutine DomainChange
	
	
	! write raydata
	subroutine WriteRayData(ray, fname)
	
		type(rayContainer), intent(in) :: ray
		character(len=*), intent(in)   :: fname
		
		open(unit=81, file=fname, action='write',  position='append') 
	    write(81,'(7(1x,e14.6))') ray%point, ray%direction, ray%power
	    close(unit=81)
	        
	end subroutine WriteRayData
	
	
	! get angle of incident
	! return cosine of incident angle and normal vector pointing inside the
	! current tetra-element
	subroutine IncidentAngle(tetra, ray, cosAngle, nsf)
	
		type(tetraElement), intent(in)      :: tetra
		type(rayContainer), intent(in)      :: ray
		real(dp), intent(out)               :: cosAngle
		real(dp), dimension(3), intent(out) :: nsf
		real(dp), dimension(3)              :: dummy
	
		! get surface normal of current face
	    call return_surfNormal(tetra, ray%faceID, nsf, dummy)
	    
	    ! cosine of angle between surface normal and incident direction
        ! must be positive (else normal is pointing outwards and must be used as negative)
        cosAngle = -dot_product(ray%direction,nsf)
        if (cosAngle < 0.0_dp) then
            write(*,*) "change direction of normal"
            cosAngle = -cosAngle
            nsf = -nsf
        end if
	
	end subroutine IncidentAngle
	
end module tracing