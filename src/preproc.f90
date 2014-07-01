module preproc

	use element

	implicit none

	type elementbins
		integer,allocatable :: bin(:)
	end type elementbins

	contains

	subroutine binelement(elno,dmax,dmin,ec,elbins)
		integer :: elno,nbins,binsize,cuPerD,sz,binno,sl(3)
		integer,allocatable :: tempbin(:)
		real(8) :: dmax(3),dmin(3),cent(3),ec(4,3)
		type(elementbins),intent(inout) :: elbins(:)

		nbins = size(elbins,1)
		cuPerD = nbins**(1.d0/3.d0)
		cent = sum(ec,1)/4.d0
		sl = ceiling(((cent-dmin)/(dmax-dmin))*cuPerD)
		binno = sl(1) + (sl(2)-1)*cuPerD + (sl(3)-1)*cuPerD**2
		sz = size(elbins(binno)%bin,1)
		allocate(tempbin(sz+1))
		tempbin(1:sz) = elbins(binno)%bin
		tempbin(sz+1) = elno
		call move_alloc(tempbin,elbins(binno)%bin)

	end subroutine binelement

	subroutine getpointbin(nbins,dmax,dmin,pt,binno)
		integer :: nbins,binno,cuPerD,sl(3)
		real(8) :: dmax(3),dmin(3),pt(3)

		cuPerD = nbins**(1/3)
		sl = ceiling(((pt-dmin)/(dmax-dmin))*cuPerD)
		binno = sl(1) + sl(2)*cuPerD + sl(3)*cuPerD**2
	end subroutine getpointbin

end module preproc
