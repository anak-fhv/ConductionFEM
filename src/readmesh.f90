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
