! Raytracing element property module module
! author: Steffen Finck
! contact: steffen.finck@fhv.at

module rt_funcs

    use rt_types
    use math_funs
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
    
!     ! return face number for (surface) elements, i.e. handle negative face IDs
!     integer function get_facenumber(face, tetra)
!     
! 	    integer            :: face
! 	    type(tetraElement) :: tetra
! 	    integer :: id
! 	    
! 	    if (face > 0) then
! 		    get_facenumber = face
! 		else if (face < 0) then
! 			if (count(face == tetra%neighbors(:,2)) /= 1) then
! 				write(*,*) "tetraelement has more than 1 or none face on the surface!"
! 				stop
! 		    else 
! 			    do id = 1,4
! 				    if (face == tetra%neighbors(id,2)) exit
! 				end do
! 				get_facenumber = id
! 		    end if
! 		else
! 			write(*,*) "face Id is 0!"
! 			stop
! 	    end if
!     
!     end function get_facenumber
    
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
    subroutine return_coords(tetra, vertices, vertIDs, p1, p2, p3)
    
	    type(tetraElement), intent(in)      :: tetra
	    integer, intent(in), dimension(3)   :: vertIDs
	    real(dp), intent(in)                :: vertices(:,:)
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
        PointInside = (all(bc > 0e0_dp) .and. (abs(sum(bc) - 1) <= 1e-14_dp))
        
        if (PointInside .eqv. .false.) then
	        write(*,*) bc
	        write(*,*) abs(sum(bc) - 1)
        end if
        
    end function
    
	! return surface normal and reference point (one of thee vertices)
	subroutine return_surfNormal(tetra, face, vertices, normal, point)
	
		type(tetraElement), intent(in)      :: tetra
		integer, intent(in)                 :: face
		real(dp), intent(in)                :: vertices(:,:)
		real(dp), dimension(3), intent(out) :: normal, point
		real(dp), dimension(3) :: p2, p3
		integer, dimension(3) :: vertIDs
! 		integer :: tmpface
		
! 		tmpface = get_facenumber(face,tetra)
		call return_facevertIds(face,vertIDs)  
        call return_coords(tetra, vertices, vertIDs, point, p2, p3)
        normal = cross(p2-point,p3-point)
        normal = normal/norm(normal)
	
	end subroutine return_surfNormal
	
end module rt_funcs
