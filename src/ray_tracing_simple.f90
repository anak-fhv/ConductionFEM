program ray_tracing_simple

	implicit none

	type elementbins
		integer,allocatable :: bin(:)
	end type elementbins

	contains

	subroutine binelement(dmax,dmin,ec,elbins)
		integer :: nbins,binsize,cuPerD,sl(3)
		real :: dmax(3),dmin(3),cent(3),ec(4,3)
		type(elementbins) :: elbins(:)

		nbins = size(elbins,1)
		cuPerD = nbins**(1/3)
		cent = sum(ec,1)/4.d0
		sl = ceiling(((cent-dmin)/(dmax-dmin))*cuPerD)
		binno = sl(1) + sl(2)*cuPerD + sl(3)*cuPerD**2
		sz = size(elbins(binno)%bin,1)
		allocate(tempbin(sz+1))
		tempbin(1:sz) = elbins(binno)%bin
		tempbin(sz+1) = elno
		call move_alloc(tempbin,elbins(binno)%bin)

	end subroutine binelement

	subroutine getpointbin(nbins,dmax,dmin,pt,binno)
		integer :: nbins,binno,sl(3)
		real(8) :: dmax(3),dmin(3),pt(3)

		cuPerD = nbins**(1/3)
		sl = ceiling(((cent-dmin)/(dmax-dmin))*cuPerD)
		binno = sl(1) + sl(2)*cuPerD + sl(3)*cuPerD**2
	end subroutine getpointbin

end program ray_tracing_simple
