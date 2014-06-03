! Raytracing tracing general math module

module rt_math
    
    use rt_parameters
    implicit none
    
    contains
    
    ! put all zeros at the end of a list
    ! list should contain only positive or zero entries
    subroutine putZerosAway(list)

        integer, intent(inout) :: list(:)
        integer                :: nmax, nzeros, value
        integer, dimension(1)  :: id
        
        nmax = size(list)
        nzeros = count(list == 0)
        if (nzeros == 0) return
 
        ! more of a debuuging feature, can be commented size
        if (minval(list) < 0) then
            write(*,*) "list should not contain negative entries!"
            write(*,*) list
            stop
        end if

        do
            id = minloc(list)                 ! first index of smallest element 
            if (nmax-nzeros == id(1)-1) exit  ! all zero's are at the end exit loop
            
            ! move zero at the end, and shift remaining list about 1 element to the front
            value = list(id(1))
            list(id(1):nmax-1) = list(id(1)+1:nmax) 
            list(nmax) = value
            
        end do
        
    end subroutine putZerosAway

    ! determine a rotation matrix given an axis of rotation (rotAxis),
    ! an angle of rotation (rotAngle), and a vector perpendicular (vecPerp) to it
    subroutine RotationMatrix(rotAxis, vecPerp, rotAngle, rotMatrix)
    
         real(dp), dimension(3), intent(in)    :: rotAxis, vecPerp
         real(dp)                              :: rotAngle
         real(dp), dimension(3,3), intent(out) :: rotMatrix
         real(dp), dimension(3,3)              :: A,M
         
         ! cartesian basis
         A = reshape([rotAxis, vecPerp, cross(rotAxis, vecPerp)],[3,3])
         ! general rotation matrix
         M = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(rotAngle), sin(rotAngle), 0.0_dp, -sin(rotAngle), cos(rotAngle)],[3,3])
         !rotation matrix for cartesian basis
         rotMatrix = matmul(matmul(A,M),transpose(A))
         
    end subroutine RotationMatrix
    
    ! cross product in 3d
    function cross(v1,v2)

        real(dp), dimension(3) :: cross
        real(dp), dimension(3), intent(in)  :: v1,v2
 
        cross(1) = v1(2)*v2(3) - v1(3)*v2(2) 
        cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
 
    end function cross
  
    ! norm 
    function norm(v)

        real(dp)             :: norm
        real(dp), intent(in) :: v(:)
  
        norm = sqrt(dot_product(v,v))

	end function norm
	
	! linear interpolation between value pairs (x1,y1) and (x2,y2)
	! output: value at given point x3
	real(dp) function linInterpol(x1,x2,y1,y2,x3)
	
		real(dp) :: x1,x2,x3,y1,y2, deltaY, deltaX
		
		deltaY = y1 - y2
        deltaX = x1 - x2
        linInterpol = deltaY/deltaX*(x3 - x1) + y1
        
        if (abs(deltaX) < 1e-13_dp) write(*,*) "numerics!"
    
    end function linInterpol
    
end module rt_math                           
