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
        real(dp)           :: t1, t2, pot, totalpower
        integer            :: io_error, write_error, k, alloc_status, i
        character(len=100) :: leaveFname, d_file1, d_file2, tmp 
        type(runstatistic) :: stats
        
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
        allocate(powerNodal(size(vertices,1)), stat=alloc_status)
        call check_alloc_error(alloc_status, "powerNodal vector")
        powerNodal = 0.0_dp

        
        ! do raytracing
        write(*,*)
        write(*,*) "Start Raytracing"
        write(*,*) "================"
        write(*,*)
        write(*,'(a10,1x,a10,1x,a14,1x,a14,1x,a14,1x,a14)') "iter", "nonzeros", "max", "mean", "rms", "var"          
        call cpu_time(t1)
        pot = 0.0
	    do k = 1,nrays
		    
		    ! creating files for raytracing output
		    if (k <= nRayPaths) then
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
	        call CreateRay(k, ray, d_file1)

	        ! assign energy to a ray and trace it
            if (k <= nRayPaths) call WriteRayData(ray, d_file2)
	        call TraceRay(ray, leaveFname, d_file2, k<=nRayPaths)
	        
	        ! catch some statistics
	        if (k .eq. nint(10**(pot))) then
	            call runStats(stats)
		        write(*,'(i10,1x,i10,1x,e14.6,1x,e14.6,1x,e14.6,1x,e14.6)') k, stats%entries, stats%maxvalue, stats%logmeannormal, stats%dist, stats%logvarnormal 
		        if (k .le. 10) pot = pot + 0.2
		        if (k .gt. 10) pot = pot + 0.1
		    end if    
	        
	        
	    end do
	    call cpu_time(t2)
    
        write(*,*)
	    write(*,*) "run time raytracing: ", t2-t1
    
    end subroutine start_tracing
    
    
    ! given an emission surface, the following subroutine performs the following steps:
    ! 1. select a face on the emission surface by roulette wheel selection
    ! 2. select a point within the triangular face based on 2 random numbers 
    ! 3. select a random direction    
    subroutine CreateRay(k, ray, fname)

        integer, intent(in)             :: k
        type(rayContainer), intent(out) :: ray
        character(len=*), intent(in)    :: fname
        integer                         :: i, j, tmp
        integer, dimension(1)           :: id 
        real(dp)                        :: theta, ratio, totalpower 
        type(tetraElement)              :: tetra
        type(emissionSurface)           :: ems
        real(dp), dimension(3,3)        :: M
        real(dp), dimension(3,4)        :: corners
        
        
        ! get emission surface and choose tetra on face
        if (RT_setup .eq. 'led') then
			
			! led setup (contains only 1 emission surface)
			ems = emSurf(1)
			
			! choose tetra and corresponding face based on random value
			! assumption: temperature at emission surface is constant and equal everywhere
			call emFace(ems, ray)
			
			! get point and direction
			tetra = tetraData(ray%tetraID)
			call emPoint(ray, tetra, ems)
			
			! ray power
	        ray%power = ems%power/nrays
			
			! wavelength through linear interpolation
	        theta = myRandom(0)
            id = minloc(abs(spectrumB(:,2) - theta))
            if (abs(spectrumB(id(1),2) - theta) .le. 1e-13_dp) then
	            ray%wavelength = spectrumB(id(1),1)
            elseif (spectrumB(id(1),2) > theta) then
		        ray%wavelength = linInterpol(spectrumB(id(1)-1,2),spectrumB(id(1),2),spectrumB(id(1)-1,1),spectrumB(id(1),1),theta)
		    else
			    ray%wavelength = linInterpol(spectrumB(id(1),2),spectrumB(id(1)+1,2),spectrumB(id(1),1),spectrumB(id(1)+1,1),theta)
			end if
			
		elseif (RT_setup .eq. 'tomo') then
			
			! emission surface selection is power ratio based 
			if (k .le. floor(emSurf(1)%power*nrays/sum(emSurf%power))) then
				ems = emSurf(1)
				if ((k==1) .or. (k==floor(emSurf(1)%power*nrays/sum(emSurf%power)))) write(*,*) k, 1, k*sum(emSurf%power)/nrays
		    elseif (k .le. floor(sum(emSurf(1:2)%power)*nrays/sum(emSurf%power))) then 
		        ems = emSurf(2)
		        if ((k==ceiling(emSurf(1)%power*nrays/sum(emSurf%power))) .or. (k==floor(sum(emSurf(1:2)%power)*nrays/sum(emSurf%power)))) write(*,*) k, 2, k*sum(emSurf%power)/nrays
			else
			    ems = emSurf(3)
			    if ((k==ceiling(sum(emSurf(1:2)%power)*nrays/sum(emSurf%power))) .or. (k==nrays)) write(*,*) k, 1, k*sum(emSurf%power)/nrays
			end if
			
			! choose tetra and corresponding face based on random value 
			! (roulette wheel selection)
			call emFace(ems, ray)
			
			! ray power
			ray%power = sum(emSurf%power)/nrays
			
			! move to other domain in case of emission from interface
			if (any(ems%name == ignoredSurfaces) .eqv. .true.) then
				tmp = ray%tetraID
				ray%tetraID = tetraData(tmp)%neighbors(ray%faceID,1)
				ray%faceID = tetraData(tmp)%neighbors(ray%faceID,2)
				
				if (tetraData(ray%tetraID)%domain .ne. 2) then
					write(*,*) "something went wrong with the domains!"
					write(*,*) tmp, tetraData(tmp)%domain
					write(*,*) ray%tetraID, tetraData(ray%tetraID)%domain
					stop
			    end if
			end if
				
			! get point and direction
			! also substract power emitted in case of emission from interface
			tetra = tetraData(ray%tetraID)
			call emPoint(ray, tetra, ems)
				
		end if
		
		Etotal = Etotal + ray%power  ! remember power emitted 
		tetra%nrays = tetra%nrays + 1
        
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
        real(dp)                          :: length, tmp
        integer                           :: cface 
        type(tetraElement)                :: tetra
        
        ! calculate length of initial ray
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
				
		    elseif ((RT_setup .eq. 'led') .and. (tetra%domain /= tetraData(tetra%neighbors(ray%faceID,1))%domain)) then		
				
				! ray is reflected into emitting surface
				! complete absorption
				call RayAbsorbing(ray, tetra, 1.0_dp)
				! write point on path
			    if (writeflag) call WriteRayData(ray, rtfname)
			    
			    ! reset ray%length and get new mean free path length
				ray%length = 0.0_dp
				length = GetPathLength()
			    
            elseif ((RT_setup .eq. 'led') .and. (ray%length >= length) ) then 
			    
			    ! setup specific for 'led'    
				! mean-free path length is reached
				
			    ! set ray point to point where interaction happens
			    ray%point = ray%point + (length - ray%length)*ray%direction
            
			    ! absorption and scattering
			    ! happens only for blue light
			    ! quantum efficiency = 0.95, i.e. 0.95 will enter absorbtion process
			    ! of these 0.05 will be completely absorbed and 0.95 will be remmitted
			    ! yellow light will be reflected completely (no self-absorbtion)
			    ! TODO: define constants as input parameters
			    if ((ray%colorchange .eqv. .false.) .and. (myRandom(0) .le. qe)) call RayAbsorbing(ray, tetra, absorptionPercentage)
			    call RayScatter(ray, tetra, leaveFname)
						
				! write point on path
			    if (writeflag) call WriteRayData(ray, rtfname)
					    
				! reset ray%length and get new mean free path length
				ray%length = 0.0_dp
				length = GetPathLength()	    
                
			elseif ((RT_setup .eq. 'tomo') .and. (tetra%domain /= tetraData(tetra%neighbors(ray%faceID,1))%domain)) then  
			
				! setup specific for 'tomo'  
				! check if domain will change (indicates reflection or absorption at interface) 		
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
            
            ! get surface normal
            call return_surfNormal(tetra, f, nsf, p1)

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
	    real(dp)                          :: t1, t2, t3, theta
	    
	    ! do some shape function magic
	    call Cartesian2Tetra(ray%point, tetra, t1, t2, t3)
	    
	    ! save absorption values
	    powerNodal(tetra%vertexIds(1)) = powerNodal(tetra%vertexIds(1)) + t1*frac*ray%power
        powerNodal(tetra%vertexIds(2)) = powerNodal(tetra%vertexIds(2)) + t2*frac*ray%power
        powerNodal(tetra%vertexIds(3)) = powerNodal(tetra%vertexIds(3)) + t3*frac*ray%power
        powerNodal(tetra%vertexIds(4)) = powerNodal(tetra%vertexIds(4)) + max((1-t1-t2-t3)-1e-13_dp,0.0_dp)*frac*ray%power
        
        ! update raypower
        ray%power = (1.0_dp-frac)*ray%power
        
        ! change to yellow light
        if (RT_setup .eq. 'led') then
	        ray%colorchange = .true.
	        call RayWavelength(ray,spectrumY)
        end if
        
    end subroutine RayAbsorbing
    
    
    ! perform scattering (so far isotropic only)
    subroutine RayScatter(ray, tetra, fname)
    
	    type(rayContainer), intent(inout) :: ray
	    type(tetraElement), intent(in)    :: tetra
	    character(len=*), intent(in)      :: fname
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
        
        ! if scattering to a boundary happens
	    if (tetra%neighbors(ray%faceID,2) < 0) then
			call BoundaryHandling(ray, tetra, fname)
	    else
	        ! output should be point on a tetraeder face with a direction inwards the tetra
	        ray%tetraID = tetra%neighbors(ray%faceID,1) ! neighbouring tetra
		    ray%faceID = tetra%neighbors(ray%faceID,2)  ! face in neighbouring tetra
        end if
        
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
		
! 		if (any(-tetra%neighbors(ray%faceID,2) == emSurf%originalID)) write(*,*) 'should leave' 
		
		if (RT_setup .eq. 'tomo') then
			if (count(-tetra%neighbors(ray%faceID,2) == emSurf%originalID) == 0) then
				! one of the bounding box boundaries
				! here specular reflection occurs
				call IncidentAngle(tetra, ray, cosAngle, nsf)
				ray%direction = ray%direction + 2*cosAngle*nsf
		    else
			    ! one of the boundary emission surfaces
			    ! here complete absortption occurs (but should not go into absorbed array)
				! call RayAbsorbing(ray, tetra, 1.0_dp)
				Eleft = Eleft - ray%power
				ray%power = 0.0_dp
		    end if
		elseif (RT_setup .eq. 'led') then
		
			! refraction indices (should come from outside)
			! refraction indices are defined by module rt_properties			
			ratio = refracIndices(tetra%domain)/refracIndices(2)
	    
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
		            
		            ! just some internal energy accounting
		            Eleft = Eleft + ray%power
		            
		            if (ray%colorchange) then
			            Eyellow = Eyellow + ray%power
			        else
			            Eblue = Eblue + ray%power
	                end if
		            ! update ray container
		            ray%power = 0.0_dp
		        end if
	        
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
		
		if (reflection .eq. 'specular') then
		
			! specular case	
			! get incident angle
			call IncidentAngle(tetra, ray, cosAngle, nsf)
        
	        ! compare with random if reflection or absorption occurs
	        if (myRandom(0) <= 1.0_dp - 1.5_dp*alpha*cosAngle) then
		        ! ray is reflected
		        ray%direction = ray%direction + 2*cosAngle*nsf  
		    else
			    ! ray is absorbed
			    ! TODO: case of partial absorption
			    call RayAbsorbing(ray, tetra, 1.0_dp)
			end if
		
		elseif (reflection .eq. 'diffuse') then
		
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
		else
			write(*,*) "UNKNOWN TYPE OF REFLECTION!"
			stop
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
	
	
	! determine wavelength in case of led setup
	subroutine RayWavelength(ray,spectrum)
	
		type(rayContainer)    :: ray
		real(dp), intent(in)  :: spectrum(:,:)
		real(dp)              :: theta
		integer, dimension(1) :: id
		
		theta = myRandom(0)
        id = minloc(abs(spectrum(:,2) - theta))
        if (abs(spectrum(id(1),2) - theta) .le. 1e-13_dp) then
	        ray%wavelength = spectrum(id(1),1)
        elseif (spectrum(id(1),2) > theta) then
		    ray%wavelength = linInterpol(spectrum(id(1)-1,2),spectrum(id(1),2),spectrum(id(1)-1,1),spectrum(id(1),1),theta)
		else
			ray%wavelength = linInterpol(spectrum(id(1),2),spectrum(id(1)+1,2),spectrum(id(1),1),spectrum(id(1)+1,1),theta)
	    end if
	
	end subroutine RayWavelength
	
	
	! calculating some runstatistics
	subroutine runStats(stats)
		type(runstatistic), intent(out) :: stats
		real(dp), dimension(:), allocatable, save :: oldValue
		
		if (allocated(oldValue) .eqv. .false.) allocate(oldValue(size(powerNodal)))
		
		stats%entries = count(abs(powerNodal).gt.1e-13_dp)
		stats%maxvalue = maxval(powerNodal)
		stats%mean = sum(powerNodal, powerNodal.gt.0)/stats%entries
		stats%logmean = sum(log10(powerNodal), powerNodal.gt.0)/stats%entries
		stats%rms = sqrt(sum(powerNodal**2,powerNodal.gt.0)/stats%entries)
		stats%var = sum((powerNodal - stats%mean)**2, powerNodal.gt.0)/(stats%entries-1)
		stats%logvar =  sum((log10(powerNodal) - stats%logmean)**2, powerNodal.gt.0)/(stats%entries-1)
		stats%logmeannormal = sum(log10(powerNodal)/log10(stats%maxvalue), powerNodal.gt.0)/stats%entries
		stats%logvarnormal =  sum((log10(powerNodal)/log10(stats%maxvalue) - stats%logmeannormal)**2, powerNodal.gt.0)/(stats%entries-1)
		stats%dist = sqrt(sum((powerNodal/stats%maxvalue-oldValue/maxval(oldValue))**2)/stats%entries)   
		oldValue = powerNodal
		
	end subroutine runStats	
	
	
	! choose emission face based on random selection
	subroutine emFace(ems,ray)
	
		type(emissionSurface), intent(in) :: ems
		type(rayContainer), intent(out)   :: ray
		real(dp)                          :: psi
		integer                           :: i
		
		! choose tetra based on random value
		! decision whether to start at the beginning or end of the list (just for speed)
	    psi = myRandom(0)
	    if (psi > 0.5) then
	        do i = size(ems%value)-1,1,-1
	            if (ems%value(i) < psi) exit
	        end do
	        ray%tetraID = ems%elemData(i+1,1)
	        ray%faceID = ems%elemData(i+1,2)
	    else           
	        do i = 1,size(ems%value)
	            if (ems%value(i) > psi) exit
	        end do
	        ray%tetraID = ems%elemData(i,1)
	        ray%faceID = ems%elemData(i,2)
	    end if
        
	end subroutine emFace
	
	
	! get point within face and choose direction
	subroutine emPoint(ray,tetra, ems)
	
		type(rayContainer), intent(inout) :: ray
		type(tetraElement), intent(in)    :: tetra
		integer, dimension(3)             :: vertIDs 
	    real(dp), dimension(3)            :: p1, p2, p3, dir1, ndir 
	    real(dp)                          :: tc1, tc2, psi, theta
	    type(emissionSurface)             :: ems
	    
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
        ndir = cross(p2-p1, p3-p1)
        ndir = ndir/norm(ndir)
        
        ! check whether normal vector points inwards or outwards
        if (dot_product(ndir, vertices(tetra%vertexIds(10-sum(vertIDs)),:) - ray%point) < 0) ndir = -ndir

        ! rotate normal vector into direction defined by psi and theta
        dir1 = (p2-p1)/norm(p2-p1)
        theta = asin(sqrt(myRandom(0)))
        psi = 2.0_dp*pi*myRandom(0) 	         
        ray%direction = sin(theta)*(cos(psi)*dir1 + sin(psi)*cross(ndir, dir1)) + cos(theta)*ndir
        
        ! substract in case of tomo data and emission from interface
        if ((any(ems%name == ignoredSurfaces) .eqv. .true.) .and. (RT_setup .eq. 'tomo')) then
	        powerNodal(tetra%vertexIds(1)) = -ray%power*tc1
	        powerNodal(tetra%vertexIds(2)) = -ray%power*tc2
	        powerNodal(tetra%vertexIds(3)) = -ray%power*(1.0_dp-tc1-tc2)
	        Etotal = Etotal - ray%power 
	    end if
	    
	end subroutine emPoint
	
end module tracing