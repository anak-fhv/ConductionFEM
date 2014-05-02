module element

	implicit none

	type elementbins
		integer,allocatable :: bin(:)
	end type elementbins

	contains

	subroutine shapefunctions(ec,v6,sp)
		real(8),dimension(4,3) :: ec
		real(8),dimension(4,4) :: sp,vm
		real(8),dimension(3,3) :: ed
		real(8) :: v6
		integer :: i,j,k

		vm(:,2:4) = ec
		vm(:,1) = 1.0d0
		v6 = abs(det4(vm))
		call invertcoords(vm)
		sp = transpose(vm);

	end subroutine shapefunctions

	subroutine invertcoords(vm)
		real(8),dimension(4,4) :: vm
		real(8),dimension(256) :: work
		integer,dimension(4) :: ipiv
		integer,parameter :: m=4,n=4,lda=4,lwork=256
		integer :: info

		call dgetrf(m,n,vm,lda,ipiv,info)
		call dgetri(m,vm,lda,ipiv,work,lwork,info)
	end subroutine invertcoords

	subroutine elementstiffness(sp,ev,elk,btdb)
		real(8),dimension(4,4) :: sp,btdb
		real(8),dimension(4,3) :: bt
		real(8),dimension(3,4) :: b
		real(8) :: ev,elk

		bt = sp(:,2:4)
		b = transpose(bt)
		btdb = elk*(matmul(bt,b))*(ev/6.0d0)

	end subroutine elementstiffness

	subroutine elementcentroid(ec,centroid)
		real(8) :: centroid(3),ec(4,3)

		centroid = sum(ec,1)/4.d0
	end subroutine elementcentroid

	subroutine addtoelementbins(elno,elcent,dnum,dlow,dhigh,elbins)
		integer :: n,nentries,elno,dnum,binnum
		integer,allocatable :: temp(:)
		real(8) :: dlow,dhigh,elcent(3)
		type(elementbins) :: elbins(:)

		n = size(elbins,1)
		binnum = ceiling(real(n)*((elcent(dnum)-dlow)/(dhigh-dlow)))
		print *, "dlow: ", dlow
		print *, "dhigh: ", dhigh
		print *, "elcent(3); ", elcent(3)
		print *, "binnum: ", binnum

		if(.not.(allocated(elbins(binnum)%bin))) then
			allocate(elbins(binnum)%bin(1))
			elbins(binnum)%bin(1) = elno
		else
			nentries = size(elbins(binnum)%bin,1)
			allocate(temp(nentries+1))
			temp(1:nentries) = elbins(binnum)%bin
			temp(nentries+1) = elno
			call move_alloc(temp,elbins(binnum)%bin)
		end if
		
	end subroutine addtoelementbins

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

	function shapefuncsquaresurfint(fcnodes) result(surfint)
		integer,dimension(3) :: fcnodes
		real(8),dimension(4,4) :: surfint
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
		real(8),dimension(3,3) :: fc
		real(8),dimension(3) :: v,w,a
		real(8) :: area

		v = fc(2,:)-fc(1,:)
		w = fc(3,:)-fc(1,:)
		a = cross_product_3(v,w)
		area = 0.5d0*norm2(a)
	end function facearea

	function cross_product_3(v,w) result(a)
		real(8),dimension(3) :: v,w,a
		integer :: i

		a(1) = v(2)*w(3)-w(2)*v(3)
		a(2) = -(v(1)*w(3)-w(1)*v(3))
		a(3) = v(1)*w(2)-w(1)*v(2)
	end function cross_product_3

	function detcreate(ec,pos1,pos2) result(vald)
		real(8),dimension(4,3) :: ec
		real(8),dimension(3,3) :: vald
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
		real(8),dimension(3,3) :: A
		real(8) :: d

		d =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
	end function det3

	function det4(A) result(d)
		real(8),dimension(4,4) :: A
		real(8) :: d

		d = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+			&
			A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+					&
			A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-					&
			A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+			&
			A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ 					&
			A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+					&
			A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+			&
			A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+					&
			A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-					&
			A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ 			&
			A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+					&
			A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

	end function det4

end module element
