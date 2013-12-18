module elementcalculations
	implicit none
	contains
	
	subroutine shapefunctions(ec,v6,sp)
		double precision,dimension(4,3) :: ec
		double precision,dimension(4,4) :: sp,vm
		double precision,dimension(3,3) :: ed
		double precision :: v6
		integer :: i,j,k

		vm(:,2:4) = ec
		vm(:,1) = 1.0d0
		v6 = abs(det4(vm))
		call invertcoords(vm)
		sp = transpose(vm);
		
	end subroutine shapefunctions

	subroutine invertcoords(vm)
		double precision,dimension(4,4) :: vm
		double precision,dimension(256) :: work
		integer,dimension(4) :: ipiv
		integer,parameter :: m=4,n=4,lda=4,lwork=256
		integer :: info

		call dgetrf(m,n,vm,lda,ipiv,info)
		call dgetri(m,vm,lda,ipiv,work,lwork,info)
	end subroutine invertcoords

	subroutine elementstiffness(sp,ev,elk,btdb)
		double precision,dimension(4,4) :: sp,btdb
		double precision,dimension(4,3) :: bt
		double precision,dimension(3,4) :: b
		double precision :: ev,elk
		
		bt = sp(:,2:4)
		b = transpose(bt)
		btdb = elk*(matmul(bt,b))*(ev/6.0d0)
		
	end subroutine elementstiffness

	subroutine bfacenodes(fcnum,fcnodes)
		integer :: fcnum
		integer,dimension(3) :: fcnodes

		if(fcnum == 1) then
			fcnodes = (/1,2,3/)
		elseif(fcnum == 2) then
			fcnodes = (/1,2,4/)
		elseif(fcnum == 3) then
			fcnodes = (/2,3,4/)
		elseif(fcnum == 4) then
			fcnodes = (/1,3,4/)
		end if
	end subroutine bfacenodes

	subroutine temperatureboundary(fcnum,Telem)
		double precision,dimension(4) :: Telem
		integer :: fcnum
	end subroutine temperatureboundary

	subroutine fluxboundary(ec,fcnodes,bv,elq)
		double precision,dimension(4,3) :: ec
		double precision,dimension(3,3) :: fc
		double precision,dimension(4) :: elq
		double precision :: fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: i,j,k,el

		elq = 0.0d0
		fc = ec(fcnodes,:)
		fcarea = facearea(fc)
		elq(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv
	end subroutine fluxboundary

	subroutine convectiveboundary(ec,fcnodes,bv,Tamb,hta,elht)
		double precision,dimension(4,3) :: ec
		double precision,dimension(4,4) :: sp,elht,surfint
		double precision,dimension(3,3) :: fc
		double precision,dimension(4) :: hta
		double precision :: Tamb,fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: fnum,i,bl

		hta = 0.0d0
		elht= 0.0d0
		fc = ec(fcnodes,:)
		do i=1,3
		end do
		fcarea = facearea(fc)
		hta(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv*Tamb
		surfint = shapefuncsquaresurfint(fcnodes)
		elht = (2.0d0*fcarea/24.0d0)*bv*surfint
	end subroutine convectiveboundary

	function shapefuncsquaresurfint(fcnodes) result(surfint)
		integer,dimension(3) :: fcnodes
		double precision,dimension(4,4) :: surfint
		integer :: i,j,n1,n2

		surfint = 0.0d0
		do i=1,3
			n1 = fcnodes(i)
			do j=1,3
				n2 = fcnodes(j)
				if(n1==n2) then
					surfint(n1,n2) = 2.0d0
				else
					surfint(n1,n2) = 1.0d0
				end if
			end do
		end do
	end function shapefuncsquaresurfint

	function facearea(fc) result(area)
		double precision,dimension(3,3) :: fc
		double precision,dimension(3) :: v,w,a
		double precision :: area
		
		v = fc(2,:)-fc(1,:)
		w = fc(3,:)-fc(1,:)
		a = cross_product_3(v,w)
		area = 0.5d0*norm2(a)
	end function facearea

	function cross_product_3(v,w) result(a)
		double precision,dimension(3) :: v,w,a
		integer :: i

		a(1) = v(2)*w(3)-w(2)*v(3)
		a(2) = -(v(1)*w(3)-w(1)*v(3))
		a(3) = v(1)*w(2)-w(1)*v(2)
	end function cross_product_3

	function detcreate(ec,pos1,pos2) result(vald)
		double precision,dimension(4,3) :: ec
		double precision,dimension(3,3) :: vald
		integer,dimension(4,3) :: order
		integer :: pos1,pos2,i,j,k

		order(1,:) = (/2,3,4/)
		order(2,:) = (/3,4,1/)
		order(3,:) = (/4,1,2/)
		order(4,:) = (/1,2,3/)

		vald = ec(order(pos1,:),:)
		if(pos2/=1) then
			vald(:,pos2-1) = 1.0d0
		end if
	end function detcreate

	function det3(A) result(d)
		double precision,dimension(3,3) :: A
		double precision :: d
		
		d =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
	end function det3

	function det4(A) result(d)
		double precision,dimension(4,4) :: A
		double precision :: d

		d =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
		A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
		A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
		A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
		A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

	end function det4

end module elementcalculations
