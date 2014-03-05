module assembly

	implicit none

	type noderow
		integer,dimension(:),allocatable :: col
		real(8),dimension(:),allocatable :: val
	end type noderow

	contains

!-------------------------------------------------------------------
! Subroutine to assemble element values into the global nodes
!-------------------------------------------------------------------
	subroutine assemble_noderows(npart,elnodes,btdb)
		type(noderow),dimension(4) :: npart
		integer,dimension(4) :: elnodes
		real(8),dimension(4,4) :: btdb
		integer,dimension(:),allocatable :: tempcol
		real(8),dimension(:),allocatable :: tempval
		integer :: i,sz

		do i=1,4
			if(.not.(allocated(npart(i)%col))) then
				allocate(npart(i)%col(4))
				allocate(npart(i)%val(4))
				npart(i)%col = elnodes
				npart(i)%val = btdb(i,:)
			else
				sz = size(npart(i)%col,1)
				allocate(tempcol(sz+4))
				allocate(tempval(sz+4))
				tempcol = (/npart(i)%col,elnodes/)
				tempval = (/npart(i)%val,btdb(i,:)/)
				call move_alloc(tempcol,npart(i)%col)
				call move_alloc(tempval,npart(i)%val)
				if(allocated(tempcol)) deallocate(tempcol)
				if(allocated(tempval)) deallocate(tempval)
			end if
		end do

	end subroutine assemble_noderows
!-------------------------------------------------------------------
! END subroutine to assemble element values
!-------------------------------------------------------------------



	subroutine collapse_noderows(node,val,col,row_ptr)
		type(noderow),dimension(:) :: node
		integer,dimension(:),allocatable :: tempcol,ind,row_ptr,col
		real(8),dimension(:),allocatable :: tempval,val
		integer :: i,numunique,numno,numval,sz

		numno = size(node,1)
		allocate(row_ptr(numno+1))
		row_ptr(1) = 1
		numval = 0
		do i=1,size(node,1)
			sz = size(node(i)%col,1)
			allocate(ind(sz))
			call indexedsort(node(i)%col,ind)
			node(i)%val = node(i)%val(ind)
			deallocate(ind)
!-------------------------------------------------------------------
! This call added to make the matrix store only the upper triangle
!-------------------------------------------------------------------
!			if(i .eq. 5) then
!				write(*,'(4i8)') node(i)%col
!			end if
!			call deletelower(i,node(i)%col,node(i)%val)
!			if(i .eq. 5) then
!				write(*,'(4i8)') node(i)%col
!			end if
!-------------------------------------------------------------------
! End call block
!-------------------------------------------------------------------
			call add_duplicates(node(i)%col,node(i)%val,numunique)
			numval = numval+numunique
			row_ptr(i+1) = row_ptr(i)+numunique
		end do

		allocate(val(numval))
		allocate(col(numval))
		do i=1,numno
			val(row_ptr(i):(row_ptr(i+1)-1)) = node(i)%val
			col(row_ptr(i):(row_ptr(i+1)-1)) = node(i)%col
! Commented for debugging
			deallocate(node(i)%col)
			deallocate(node(i)%val)
! Commented for debugging
		end do

	end subroutine collapse_noderows

!-------------------------------------------------------------------
! Smaller subroutines/helper routines
!-------------------------------------------------------------------
	subroutine indexedsort(a,b)
		integer,dimension(:) :: a
		integer,dimension(size(a,1)) :: b
		integer :: i,j,n,temp

		n = size(a,1)
		b = (/1:n/)
		do i=1,n-1
			do j=i+1,n
				if(a(i) > a(j)) then
					temp = a(j)
					a(j) = a(i)
					a(i) = temp
					temp = b(i)
					b(i) = b(j)
					b(j) = temp
				end if
			end do
		end do
		temp = 0
	end subroutine indexedsort

	subroutine add_duplicates(nodecol,nodeval,numunique)
		integer,dimension(:),allocatable :: nodecol,ind,tempcol
		integer,dimension(:),allocatable :: starts
		real(8),dimension(:),allocatable :: nodeval,tempval
		integer :: i,j,k,numunique

		if(size(nodecol,1)==1) then
			numunique = 1
			return
		end if

		call find_duplicates(nodecol,starts,numunique)
		allocate(tempcol(numunique))
		allocate(tempval(numunique))
		tempcol = nodecol(starts)

		do i=1,numunique
			tempval(i) = sum(nodeval(starts(i):starts(i+1)-1))
		end do
		call move_alloc(tempcol,nodecol)
		call move_alloc(tempval,nodeval)
	end subroutine add_duplicates

	subroutine find_duplicates(a,starts,numunique)
		integer,dimension(:) :: a
		integer,dimension(:),allocatable :: starts,temp
		integer,dimension(1) :: ml
		integer :: i,j,ind,numunique

		allocate(starts(size(a,1)))
		starts = 0
		ind=1
		starts(1) = 1
		i = 2
		do while(a(ind).lt.maxval(a))
			ml = minloc(a,a.gt.a(ind))
			starts(i) = ml(1)
			ind = ml(1)
			i = i+1
		end do
		allocate(temp(i))
		temp = (/starts(1:i-1),size(a,1)+1/)
		call move_alloc(temp,starts)
		numunique = i-1
	end subroutine find_duplicates

	subroutine deletelower(n,nodecol,nodeval)
		integer :: n,i,sz,ct
		integer,allocatable :: nodecol(:),tempcol(:)
		real(8),allocatable :: nodeval(:),tempval(:)

		sz = size(nodecol,1)
		if(sz .eq. 1) then
			return
		end if

		ct = 0
		do i=1,sz
			if(nodecol(i) .lt. n) then
				ct = ct + 1
			end if
		end do

		allocate(tempcol(sz-ct))
		allocate(tempval(sz-ct))

		tempcol = nodecol(ct+1:sz)
		tempval = nodeval(ct+1:sz)

		call move_alloc(tempcol,nodecol)
		call move_alloc(tempval,nodeval)

	end subroutine deletelower
!-------------------------------------------------------------------
! END smaller subroutines/helper routines
!-------------------------------------------------------------------

end module assembly
