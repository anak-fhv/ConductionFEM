! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 21.02.2014
! version: 0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
module array
    use rt_types
    
    
    
end module

program raytracing

    use math_funs
    use pre_process_data
    use rt_funcs
    use array
   
    implicit none
    character(len=100)                               :: file_name
    character(len=100), dimension(1)                 :: emSurfNames
    real(dp), dimension(3)                           :: rayOrigin, rayDir
    integer                                          :: npart, i, n, k, alloc_status
    real                                             :: t1, t2
    integer, dimension(:), allocatable               :: counter
    real(dp)                                         :: psi, maxratio, minratio
    integer                                          :: kmax = 10000000
    type(tetraElement), dimension(:), allocatable    :: tetras
    type(emissionSurface), dimension(:), allocatable :: emSurf
    real(dp), dimension(:,:), allocatable            :: vertexData
    
    ! user input
    file_name = "a.msh"
    npart = 1
    emSurfNames = ['xLow'] ! right now this requires to state the correct dimension in above declaration       
    psi = myRandom(2903)
    
    
    ! perform pre-processing
    call read_mesh_data(file_name, emSurfNames, npart, tetras, vertexData, emSurf)
    
    call cpu_time(t1)
    do i = 1,2
        call CreateRay(tetras, vertexData, emSurf(1), rayOrigin, rayDir)
    end do
    call cpu_time(t2)
    
    write(*,*) "time: ", t2-t1

!        
!     ! test whether new point is on face of new tetra-element
!     face = tetras(nTetra)%neighbors(newFace,2)
!     nTetra = tetras(nTetra)%neighbors(newFace,1)
!   
!     call return_facevertIds(face,vertIds)  
!     p1 = vertices(:,tetras(nTetra)%vertexIds(vertIds(1)))
!     p2 = vertices(:,tetras(nTetra)%vertexIds(vertIds(2)))
!     p3 = vertices(:,tetras(nTetra)%vertexIds(vertIds(3)))
!   
!     na = cross(p2-p1,p3-p1)
!     nx = cross(p3-p2,waypoint-p2)
!     ny = cross(p1-p3,waypoint-p3)
!     nz = cross(p2-p1,waypoint-p1)
!     bc(1) = dot_product(na,nx)
!     bc(2) = dot_product(na,ny)
!     bc(3) = dot_product(na,nz)
!   
!     write(*,*) bc

 contains 
 
     subroutine CreateRay(tetraData, vertices, ems, rayOrigin, rayDir)

    ! procedure for creating starting point for ray on emission surface ems
    ! procedure chooses a random tetraeder on the surface and the creates a
    ! random point with the face on the surface
        use math_funs
        implicit none
    
        type(tetraElement)                  :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(emissionSurface), intent(in)   :: ems
        real(dp), dimension(3), intent(out) :: rayOrigin, rayDir
        integer                             :: face, tetra, i
        integer, dimension(3)               :: vertIDs 
        real(dp)                            :: psi, theta, b, c, d
        real(dp), dimension(3)              :: p1, p2, p3, dir1, dir2, dir21, ds1, ndir
        real(dp), dimension(3,3)            :: A, M  
    
        ! get random tetraeder on the emission surface
        psi = myRandom(0)
        ! decision whether to start at the beginning or end of the list just for speed
        if (psi > 0.5) then
            do i = size(ems%area)-1,1,-1
                if (ems%area(i) < psi) exit
            end do
            tetra = ems%elemData(i+1,1)
            face = ems%elemData(i+1,2)
        else           
            do i = 1,size(ems%area)
                if (ems%area(i) > psi) exit
            end do
            tetra = ems%elemData(i,1)
            face = ems%elemData(i,2)
        end if
    
!         write(*,*) size(tetraData)
!         write(*,*) size(vertices)
!         write(*,*) face
!         write(*,*)
    
        ! get points of the face on the surface
        call return_facevertIds(face,vertIDs)  
        p1 = vertices(tetraData(tetra)%vertexIds(vertIDs(1)),:)
        p2 = vertices(tetraData(tetra)%vertexIds(vertIDs(2)),:)
        p3 = vertices(tetraData(tetra)%vertexIds(vertIDs(3)),:)
     
        ! find longest side of triangle
        if ((norm(p2-p1) >= norm(p3-p1)) .and. (norm(p2-p1) >= norm(p3-p2))) then
            dir1 = p2-p1
            dir2 = p3-p1
            rayOrigin = p1     
!             write(*,*) "p1"
!             write(*,*)   
        else if ((norm(p1-p3) >= norm(p2-p1)) .and. (norm(p1-p3) >= norm(p3-p2))) then
            dir1 = p1-p3
            dir2 = p2-p3 
            rayOrigin = p3
!             write(*,*) "p3"
!             write(*,*)  
        else
            dir1 = p3-p2
            dir2 = p1-p2
            rayOrigin = p2
!             write(*,*) "p2"
!             write(*,*)
        end if
        
        ! projection of dir2 onto dir1
        dir21 = dot_product(dir1,dir2)/dot_product(dir1,dir1) * dir1
    
!         write(*,*) b
!         write(*,*) c
!         write(*,*) c/b
!         write(*,*)
!         write(*,*) rayOrigin
    
        ! choose random numbers based on triangle distribution
        psi = myRandom(0)
        theta = myRandom(0)
!         write(*,*) psi

        b = norm(dir1)     ! length of dir1
        c = norm(dir21)    ! length of projection dir21
        ds1 = dir2 - dir21 ! perpendicular to d1 and going through vertex of triangle which is not on d1

        if (psi <= c/b) then
            rayOrigin = rayOrigin + sqrt(psi*b*c)*dir1/b                       ! random number along longest side (triangle distribution)
            d = sqrt(psi*b*c)*tan(acos(dot_product(dir1,dir2)/(b*norm(dir2)))) ! length of perpendicular vector at new rayOrigin
            rayOrigin = rayOrigin + theta*d*ds1/norm(ds1)                        ! random number along perpendicular direction (uniform distribution)     
        else
            rayOrigin = rayOrigin + (b - sqrt(b*(b-c)*(1-psi)))*dir1/b                               ! random number along longest side (triangle distribution)
            d = sqrt(b*(b-c)*(1-psi)) *tan(acos(dot_product(-dir1,(dir2-dir1))/(b*norm(dir2-dir1)))) ! length of perpendicular vector at new rayOrigin
            rayOrigin = rayOrigin + theta*d*ds1/norm(ds1)                                              ! random number along perpendicular direction (uniform distribution)        
        end if
        
        ! test
        !write(*,*) PointInside(p1,p2,p3,rayOrigin)
    
        ! choose ray direction
        ndir = cross(dir1,dir2) ! normal vector of triangle
        ndir = ndir/norm(ndir)
        ! check whether normal vector points inwards or outwards
        do i = 1,4
            if (any(i == vertIDs) .eqv. .false.) exit
        end do        
        if (dot_product(ndir, vertices(tetraData(tetra)%vertexIds(i),:) - rayOrigin) < 0) ndir = -ndir

        ! first rotate about about a vector in the plane
        dir1 = dir1/b
        A = reshape([dir1, ndir, cross(dir1,ndir)],[3,3])
        theta = asin(sqrt(myRandom(0)))
        M = matmul(A,reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(theta), sin(theta), 0.0_dp, -sin(theta), cos(theta)],[3,3]))
        M = matmul(M,transpose(A))
        rayDir = matmul(M,ndir)
!         write(*,*)
!         write(*,*) rayDir
!         write(*,*)
!         write(*,*) A(1,:)
!         write(*,*) A(2,:)
!         write(*,*) A(3,:)
        
        ! second rotation with normal vector as rotation axis
        psi = 2*pi*myRandom(0)
        A = reshape([ndir, dir1, cross(ndir,dir1)],[3,3])
        M = matmul(A,reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(psi), sin(psi), 0.0_dp, -sin(psi), cos(psi)],[3,3]))
        M = matmul(M,transpose(A))
        rayDir = matmul(M,rayDir)
        
!         write(*,*)
!         write(*,'(a12,3(1x,3e14.6))') "origin: ", rayOrigin
!         write(*,'(a12,3(1x,3e14.6))') "direction: ", rayDir
!         write(*,*)
!         
!         write(*,'(a12,3(1x,3e14.6))') "dir1: ", dir1
!         write(*,'(a12,3(1x,3e14.6))') "dir2: ", dir2
!         write(*,'(a12,3(1x,3e14.6))') "ndir: ", ndir
!         write(*,'(a12,3(1x,3e14.6))') "theta: ", theta
!         write(*,'(a12,3(1x,3e14.6))') "psi: ", psi
        
    end subroutine CreateRay
    
    function PointInside(p1,p2,p3,tp)
        
        real(dp), dimension(3), intent(in) :: p1,p2,p3,tp
        logical                            :: PointInside
        real(dp), dimension(3)             :: na,nx,ny,nz,bc
        
        na = cross(p2-p1,p3-p1)
        nx = cross(p3-p2,tp-p2)
        ny = cross(p1-p3,tp-p3)
        nz = cross(p2-p1,tp-p1)
        bc(1) = dot_product(na,nx)
        bc(2) = dot_product(na,ny)
        bc(3) = dot_product(na,nz)
        PointInside = (all(bc > 0e0_dp) .and. (sum(bc) <= 1))
    
    end function
  
  
end program raytracing 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutines for raytracing 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

