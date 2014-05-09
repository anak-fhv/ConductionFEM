! Raytracing elements operation module

module rt_elements

    use rt_math 
    implicit none
    
    contains
    
    ! The function return_facenumber returns the ID of the face of a tetra
    ! given that a mask exist which indicates which vertices define the face.        
    integer function return_facenumber(mask)
            
        logical, dimension(4) :: mask
        
        if (all(mask == [.true., .true., .true., .false.])) then
            return_facenumber = 1
        else if (all(mask == [.true., .true., .false., .true.])) then
            return_facenumber = 2
        else if (all(mask == [.false., .true., .true., .true.])) then
            return_facenumber = 3
        else if(all(mask == [.true., .false., .true., .true.])) then
            return_facenumber = 4
        else
            write(*,*) "unknown mask. Can not determine correct face number!"
            stop
        end if
            
    end function return_facenumber
    
    ! returns the indices of tetraeder points (between 1 and 4) for a given face 
    subroutine return_facevertIds(face, vertIds)
    
        integer, intent(in)                :: face
        integer, dimension(3), intent(out) :: vertIds
        
        select case(face)
            case(1)
                vertIds = [1,2,3]
            case(2)
                vertIds = [1,2,4]
            case(3)
                vertIds = [2,3,4]
            case(4)
                vertIds = [1,3,4]
            case default
                write(*,*) "Unknown face number!", face
                stop
       end select
       
    end subroutine return_facevertIds    
    
    ! returns the points (in x,y,z coordinates) of a tetraeder face
    subroutine return_coords(tetra, vertIDs, p1, p2, p3)
    
	    type(tetraElement), intent(in)      :: tetra
	    integer, intent(in), dimension(3)   :: vertIDs
	    real(dp), dimension(3), intent(out) :: p1, p2, p3
	    
	    p1 = vertices(tetra%vertexIds(vertIDs(1)),:)
        p2 = vertices(tetra%vertexIds(vertIDs(2)),:)
        p3 = vertices(tetra%vertexIds(vertIDs(3)),:)
        
    end subroutine return_coords
    
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
        PointInside = (all(bc > 0e0_dp) .and. (abs(sum(bc) - 1) <= 1e-13_dp))
        
        if (PointInside .eqv. .false.) then
	        write(*,*) bc
	        write(*,*) abs(sum(bc) - 1)
        end if
        
    end function
    
	! return surface normal and reference point (one of thee vertices)
	subroutine return_surfNormal(tetra, face, normal, point)
	
		type(tetraElement), intent(in)      :: tetra
		integer, intent(in)                 :: face
		real(dp), dimension(3), intent(out) :: normal, point
		real(dp), dimension(3) :: p2, p3
		integer, dimension(3) :: vertIDs
! 		integer :: tmpface
		
! 		tmpface = get_facenumber(face,tetra)
		call return_facevertIds(face,vertIDs)  
        call return_coords(tetra, vertIDs, point, p2, p3)
        normal = cross(p2-point,p3-point)
        normal = normal/norm(normal)
	
	end subroutine return_surfNormal
	
	! transform cartesian coordiantes into tetrahedral coordinates
	subroutine Cartesian2Tetra(point, tetra, t1, t2, t3)
	
		real(dp), dimension(3), intent(in) :: point
		type(tetraElement), intent(in)     :: tetra
		real(dp), intent(out)              :: t1,t2,t3
		real(dp), dimension(3,4)           :: corners
		integer                            :: i, j
		real(dp), dimension(3,3)           :: M, Minv
		real(dp)                           :: det
		
		! get tetraeder vertex points
		forall(i = 1:4) corners(:,i) = vertices(tetra%vertexIds(i),:)
		
		! assemble matrix describing transformation from tetra-coords to cartesian-coords
	    do  i = 1,3
		    do j = 1,3
			    M(j,i) = corners(j,i) - corners(j,4) 
			end do
	    end do
	    
	    ! determine determinante
	    det = dot_product(M(:,1), cross(M(:,2), M(:,3))) 
	     
	    ! get inverse
	    Minv(1,1) = M(2,2)*M(3,3) - M(2,3)*M(3,2)
	    Minv(2,1) = M(2,3)*M(3,1) - M(2,1)*M(3,3)
	    Minv(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
	    Minv(1,2) = M(1,3)*M(3,2) - M(1,2)*M(3,3)
	    Minv(2,2) = M(1,1)*M(3,3) - M(1,3)*M(3,1)
	    Minv(3,2) = M(1,2)*M(3,1) - M(1,1)*M(3,2)
	    Minv(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
	    Minv(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
	    Minv(3,3) = M(1,1)*M(2,2) - M(1,2)*M(2,1)
	    Minv = Minv/det
	    
	    ! get tetrahedral coordinates
	    t1 = dot_product(Minv(1,:), (point - corners(:,4)))
	    t2 = dot_product(Minv(2,:), (point - corners(:,4)))
	    t3 = dot_product(Minv(3,:), (point - corners(:,4))) 
	    
	    ! numerics
	    if (abs(t1) <= 1e-13_dp) t1 = 0
	    if (abs(t2) <= 1e-13_dp) t2 = 0
	    if (abs(t3) <= 1e-13_dp) t3 = 0
	    
	    ! sanity checks
	    if ( (t1+t2+t3 > 1.0_dp + 1e-13_dp) .or. (min(t1,t2,t3) < 0) )then
		      write(*,*) "point not in tetra"
		      write(*,*) t1, t2, t3
		      write(*,*) point
		      write(*,*) corners(:,1)
	          write(*,*) corners(:,2)
	          write(*,*) corners(:,3)
	          write(*,*) corners(:,4)
	          stop
	    end if
	    
	end subroutine Cartesian2Tetra
	
end module rt_elements
