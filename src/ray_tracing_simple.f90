program ray_tracing_simple

	implicit none

!--------------------------------------------------------------------
!	Program contents

!	End program contents
!--------------------------------------------------------------------

	contains
!--------------------------------------------------------------------
!	Ray tracing main module
!--------------------------------------------------------------------
module raytracing
	use readmesh
	use preproc
	implicit none
	contains

	subroutine startray(stpt,stdir)
		
	end subroutine startray
		
	subroutine getpointbin(nbins,dmax,dmin,pt,binno)
		integer :: nbins,binno,sl(3)
		real(8) :: dmax(3),dmin(3),pt(3)

		cuPerD = nbins**(1/3)
		sl = ceiling(((cent-dmin)/(dmax-dmin))*cuPerD)
		binno = sl(1) + sl(2)*cuPerD + sl(3)*cuPerD**2
	end subroutine getpointbin

	subroutine locatepointelem(pointbin,noVerts,connTab,elem)
		integer :: i,j,k,n,elem,binelems(:)
		type(elementbins) :: pointbin,adjbins(26)

		binelems = pointbin%bin
		n = size(binelems,1)

		elem = 0	! Initialisation needed for non-bin elements

		do i=1,n
			elem = binelems(i)
			ec = noVerts(connTab(elem),:)
			call ptintettest(pt,ec,ptin)
			if(ptin) return
		end do
	end subroutine locatepointelem

	subroutine ptintettest(pt,ec,ptin)
		integer :: i
		real(8) :: sig,det,pt(3),ec(4,3),vm(4,4)
		logical :: ptin = .false.

		do i=1,5
			vm(:,1:3) = ec
			if(i.gt.1) vm(i-1,1:3) = pt
			vm(:,4) = 1
			det = det4(vm)
			sig = sig*det
		end do

		if(sig.ge.0.d0) ptin = .true.
	end subroutine ptintettest

	subroutine getbinadjacency(binno,nbins,adj)
		integer :: binno,nbins,cpd,init(2),relc(2,8),adj(3,3,3)

		cpd = nbins**(1.d0/3.d0)
		adj = 0
		init = (/1,cpd/)
		relc(1,:) = (/1,1,0,-1,-1,-1,0,1/)
		relc(2,:) = (/0,1,1,1,0,-1,-1,-1/)
		do j=1,8
			adj(2+relc(1,j),2+relc(2,j),i) = binno+i*cpd**2 + dot_product(init,relc(:,j))
		end do

		if(mod)
	end subroutine getbinadjacency

	function det4(A) result(d)
		real(8) :: d,A(4,4)

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

end module raytracing
!--------------------------------------------------------------------
!	End ray tracing module
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Read mesh module
!--------------------------------------------------------------------
module readmesh

	use preproc
    implicit none

	type surfaceData
		character(len=16) :: sfName
		integer :: sfNumber
		integer,allocatable :: elNum(:),fcNum(:)
	end type surfaceData

    contains

    subroutine openmeshfile(fno,fname)
        integer :: fno
        character(*) :: fname

        open(unit=fno,file=fname,form='formatted',status='old')
    end subroutine openmeshfile

    subroutine closemeshfile(fno)
        integer :: fno

        close(fno)
    end subroutine closemeshfile

    subroutine readmeshdetails(fno,mDets)
        integer :: fno,mDets(7)

        read(fno, *)
        read(fno, *)
        read(fno, *) mDets
    end subroutine readmeshdetails

    subroutine readmeshvertices(fno,mDets,noVerts)
        integer :: i,fno,mDets(7)
        real(8),allocatable :: noVerts(:,:)

        allocate(noVerts(mDets(1), 3))
        do i=1,mDets(1)
            read(fno, *) noVerts(i,:)
        end do
    end subroutine readmeshvertices

    subroutine readmeshconnectivity(fno,mDets,noVerts,connTab,elbins)
        integer :: fno,i,nbins,mDets(7)
        integer,allocatable :: connTab(:,:)
		real(8) :: dmax(3),dmin(3),ec(4,3),noVerts(:,:)
		type(elementbins),allocatable :: elbins(:)

		nbins = (size(noVerts,1))**(1.d0/3.d0)
		nbins = nbins**3
		do i=1,3
			dmax(i) = maxval(noVerts(:,i))
			dmin(i) = minval(noVerts(:,i))
		end do
        allocate(connTab(mDets(2),4))
		allocate(elbins(nbins))
        do i=1,mDets(2)
            read(fno, *) connTab(i,:)
			ec = noVerts(connTab(i,:),:)
			call binelement(i,dmax,dmin,ec,elbins)
        end do
    end subroutine readmeshconnectivity

    subroutine readmeshdomains(fno,mDets,doElems)
        integer :: fno,doSize,cSize,i,j,k,mDets(7),temp(10)
        integer,allocatable :: doElems(:)

        allocate(doElems(mDets(2)))
        do i=1,mDets(6)
            read(fno,'(i8)') doSize

            if(mod(doSize,10) == 0) then
                cSize = doSize/10
            else
                cSize = 1 + (doSize/10)
            end if
            do j=1,cSize
                read(fno,'(10(1x,i8))') temp
                do k=1,10
                    if(temp(k) /= 0) then
                        doElems(temp(k)) = i
                    end if
                end do
            end do
        end do
    end subroutine readmeshdomains

    subroutine readmeshsurfaces(fno,mDets,sfaces,surfNames,surfaces)
        integer :: fno,nSurf,nElems,numFc,cSize,ct,i,j,k,	&
		mDets(7),temp(10)
        integer, allocatable :: sfaces(:,:)
        character(len=16), allocatable :: surfNames(:)
		type(surfaceData),allocatable :: surfaces(:)

		nElems = mDets(2)
		nSurf = mDets(7)
        allocate(surfNames(nSurf))
        allocate(sfaces(nElems,4))
		allocate(surfaces(nSurf))
        sfaces = 0
        do i=1,mDets(7)
            read(fno,'(i8,1x,a)') numFc, surfNames(i)
			surfaces(i)%sfName = surfNames(i)
			allocate(surfaces(i)%elNum(numFc))
			allocate(surfaces(i)%fcNum(numFc))			
            if(mod(numFc,5) == 0) then
                cSize = numFc/5
            else
                cSize = 1 + numFc/5
            end if
			ct = 0
            do j=1,cSize
                read(fno, '(5(2x,i8,1x,i1))') temp
                do k=1,5
                    if(temp(2*k-1) /= 0) then
                        sfaces(temp(2*k-1),temp(2*k)) = i
						ct = ct+1
                        surfaces(i)%elNum(ct) = temp(2*k-1)
						surfaces(i)%fcNum(ct) = temp(2*k)
                    end if
                end do
            end do
        end do
    end subroutine readmeshsurfaces

	subroutine readSurfaces(fno,mDets,surfaces)
		integer :: i,j,k,numSurf,numFc,ct,cs,fno,mDets(7),temp(10)
		type(surfaceData),allocatable :: surfaces(:)

		numSurf = mDets(7)
		allocate(surfaces(mDets(7)))
		do i=1,numSurf
			surfaces%sfNumber = i
            read(fno,'(i8,1x,a)') numFc, surfaces(i)%sfName
			allocate(surfaces(i)%elNum(numFc))
			allocate(surfaces(i)%fcNum(numFc))
			if(mod(numFc,5) == 0) then
                cs = numFc/5
            else
                cs = 1 + numFc/5
            end if
			ct = 0
            do j=1,cs
                read(fno,'(5(2x,i8,1x,i1))') temp
                do k=1,5
                    if(temp(2*k-1) /= 0) then
						ct = ct+1
                        surfaces(i)%elNum(ct) = temp(2*k-1)
						surfaces(i)%fcNum(ct) = temp(2*k)
                    end if
                end do
            end do
		end do
	end subroutine readSurfaces

end module readmesh
!--------------------------------------------------------------------
!	End read mesh module
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Preprocessor module
!--------------------------------------------------------------------
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

end module preproc
!--------------------------------------------------------------------
!	End preprocessor module
!--------------------------------------------------------------------

end program ray_tracing_simple
