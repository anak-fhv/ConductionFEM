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
        character(len=100) :: resFname
        
        ! initialize random generator
        dummy = myRandom(2903)
        
        ! create file for output
        resFname = objFolder//trim(file_name)//".res"
        open(unit=81, file=resFname, action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating file for results 1",81)
        write(81,'(1x,a8,4(1x,a14))',iostat=write_error) "tetraID", "x", "y", "z", "absorbed"
        call check_io_error(write_error,"writing header results file",81)
        close(unit=81, iostat=io_error)
        call check_io_error(io_error,"closing file for results 1",81)
        
        ! additional debbunging files
        open(unit=82, file=objFolder//trim(file_name)//"-initrayloc.res", action='write', status='replace', iostat=io_error)       
        call check_io_error(io_error,"creating debug-file 1",82)
        write(82,'(3(a14))',iostat=write_error) "x", "y", "z"
        call check_io_error(write_error,"writing debug-file 1 header",82)
        close(unit=82, iostat=io_error)
        call check_io_error(io_error,"closing debug-file 1",82)
        
        ! do raytracing            
        ! in case of several emission surface a way of selecting a surface must be chosen     
	    call cpu_time(t1)
	    do k = 1,nrays
	        call CreateRay(tetraData, vertices, ems(1), ray, file_name)
	        call TraceRay(tetraData, vertices, ray, resFname)
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

        ! first rotate about about a vector in the plane
        dir1 = dir1/b
        call RotationMatrix(dir1, ndir, asin(sqrt(myRandom(0))), M1)
    
        ! second rotation with normal vector as rotation axis
        call RotationMatrix(ndir, dir1, 2.0_dp*pi*myRandom(0), M2)
        
        ! get final direction
        ray%direction = matmul(M2,matmul(M1,ndir))
        
        ! just for checking (could be commented)
        open(unit=83, file=objFolder//trim(file_name)//"-initrayloc.res", action='write', position='append')  
        write(83,'(3e14.6)') ray%point
        close(unit=83)         
        
    end subroutine CreateRay
    
    ! main raytracing routine
    subroutine TraceRay(tetraData, vertices, ray, resFname)
    
        type(tetraElement), intent(inout)   :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        character(len=*), intent(in)        :: resFname
        real(dp) :: kappa, sigma, beta, omega, length
        integer, dimension(3) :: vertIDs 
        real(dp), dimension(3) :: p1, p2, p3, endpoint
        logical :: test
        integer :: counter=0, io_error, write_error

        ! properties of the medium
        kappa = 10.0_dp
        sigma = 0.75_dp
        
        ! calculate length of initial ray
        beta = kappa + sigma        
        omega = sigma/beta
        
        ! trace ray until it leaves the enclosure or the energy is below a treshold  
        do 
	        length = 1.0_dp/beta*log(1/myRandom(0))
	        endpoint = ray%point + length*ray%direction
!             write(*,*) length
!             write(*,*) endpoint
            
            ! trace path 
            do while (ray%length < length)       
		        call FindNextTetra(vertices,tetraData,ray)
		        ! point on boundary surface
	            if (ray%faceID == -1) then
! 		            write(*,*) "ray is escaped"
			        ! do boundary handling
			        return
			    end if
            end do
            
            ! absorb power and scatter ray
            if (ray%power > 0.0001_dp) then
	            call RayAbsorbing(tetraData(ray%tetraID),ray%tetraID,endpoint,(1.0_dp-omega)*ray%power,resFname)
! 	        	call RayScatter
	            ray%power = ray%power*omega
		        ray%length = 0.0_dp
            else
	            call RayAbsorbing(tetraData(ray%tetraID),ray%tetraID,ray%point,ray%power,resFname)
                return
            end if
        end do 
        
    end subroutine TraceRay
    
    ! find next tetra Element along ray direction
    subroutine FindNextTetra(vertices,tetraData,ray)
        
        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        real(dp) :: alpha, temp
        integer  :: f, newFace
        real(dp), dimension(3) :: rp, v1, v2, nsf
        integer, dimension(3)  :: vertIDs
        
        alpha = 100.0_dp  ! some large initial value for alpha
        do f = 1,4
        
            if (f == ray%faceID) cycle  ! is face where ray is emitted
            
            ! 1 vertex and 2 vectors of current face
            call return_facevertIds(f,vertIDs)
            call return_coords(tetraData(ray%tetraID), vertices, vertIDs, rp, v1, v2)
            v1 = v1 - rp ! make a vector with origin rp  
            v2 = v2 - rp ! make a vector with origin rp
            
            ! normal vector of current face  
            nsf = cross(v1,v2)
            nsf = nsf/norm(nsf) 

            ! calculate length until intersection
            temp = (dot_product(nsf,rp) - dot_product(nsf,ray%point))/dot_product(nsf, ray%direction)
            if (temp < 0) cycle ! length can not be negative
            
            ! check if temp is the shortest length so far
            ! if true, current face is the face the ray will intersect
            if (temp < alpha) then
                alpha = temp
                newFace = f
            end if
            
        end do
        
        ! update ray container
        ray%length = ray%length+alpha
        ray%point = ray%point + alpha*ray%direction
        ray%faceID = tetraData(ray%tetraID)%neighbors(newFace,2)  
        ray%tetraID = tetraData(ray%tetraID)%neighbors(newFace,1)     
    
    end subroutine
        
    subroutine RayAbsorbing(tetra, id, point, power, resFname)
    
	    type(tetraElement), intent(inout)  :: tetra
	    integer, intent(in)                :: id
	    real(dp), dimension(3), intent(in) :: point
	    real(dp), intent(in)               :: power
	    character(len=*), intent(in)       :: resFname
	    
	    tetra%absorbed = tetra%absorbed + power
	    
	    open(unit=84, file=resFname, action='write', position='append')  
	    write(84,'(1x,i8,1x,3(e14.6,1x),e14.6)') id, point, tetra%absorbed
        close(unit=84)
        
    end subroutine RayAbsorbing
    
    ! calculate refraction with Snell's law
    subroutine SnellsLaw(incident, n1, n2, nsf, reflect, refract)
        
        real(dp), intent(in)                  :: n1, n2   ! refractive indices
        real(dp), dimension(3), intent(in)    :: incident ! incident direction 
        real(dp), dimension(3), intent(inout) :: nsf      ! normal vector
        real(dp), dimension(3), intent(out)   :: refract, reflect ! refraction and reflection direction
        real(dp) :: cosAngle, ratio
        
        cosAngle = -dot_product(incident,nsf)
        if (cosAngle < 0) then
            cosAngle = -cosAngle
            nsf = -nsf
        end if
        
        reflect = incident + 2*cosAngle*nsf
        
        ratio = n1/n2
        refract = ratio*incident + ratio*cosAngle*nsf - sqrt(1 - ratio**2*(1-cosAngle**2))*nsf
   end subroutine SnellsLaw
        
end module tracing