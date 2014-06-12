module readmesh

    implicit none

	type surfaceData
		character(len=16) :: sfName
		integer :: sfNumber
		integer,allocatable :: elNum(:),fcNum(:)
	end type surfaceData

    contains

    subroutine openmeshfile(unum,filename)
        integer :: unum
        character(len = *) :: filename

        open(unit=unum,file=filename,form='formatted',status='old')
    end subroutine openmeshfile

    subroutine closemeshfile(unum)
        integer :: unum

        close(unum)
    end subroutine closemeshfile

    subroutine readmeshdetails(unum, meshdetails)
        integer :: unum
        integer, dimension(7) :: meshdetails

        read(unum, *)
        read(unum, *)
        read(unum, *) meshdetails
    end subroutine readmeshdetails

    subroutine readmeshvertices(unum, meshdetails, vertices)
        integer :: unum
        integer, dimension(7) :: meshdetails
        integer :: i
        real(8), dimension(:,:), allocatable :: vertices

        allocate(vertices(meshdetails(1), 3))
        do i=1,meshdetails(1)
            read(unum, *) vertices(i,:)
        end do
    end subroutine readmeshvertices

    subroutine readmeshconnectivity(unum,meshdetails,connectivity)
        integer :: unum,i
        integer, dimension(7) :: meshdetails
        integer, dimension(:,:), allocatable :: connectivity

        allocate(connectivity(meshdetails(2), 4))
        do i=1,meshdetails(2)
            read(unum, *) connectivity(i,:)
        end do
    end subroutine readmeshconnectivity

    subroutine readmeshdomains(unum,meshdetails,domainelements)
        integer :: unum,domainsize,countersize,i,j,k
        integer, dimension(7) :: meshdetails
        integer, dimension(10) :: tempmat
        integer, dimension(:), allocatable :: domainelements

        allocate(domainelements(meshdetails(2)))
        do i=1,meshdetails(6)
            read(unum,'(i8)') domainsize

            if(mod(domainsize,10) == 0) then
                countersize = domainsize/10
            else
                countersize = 1 + (domainsize/10)
            end if
            do j=1,countersize
                read(unum,'(10(1x,i8))') tempmat
                do k=1,10
                    if(tempmat(k) /= 0) then
                        domainelements(tempmat(k)) = i
                    end if
                end do
            end do
        end do
    end subroutine readmeshdomains

    subroutine readmeshsurfaces(unum,meshdetails,sfaces,surfacenames)
        integer :: unum,surfnum,countersize,i,j,k,currentsize,		&
		newsize, oldsize
        integer, dimension(7) :: meshdetails
        integer, dimension(10) :: tempmat
        integer, dimension(:,:), allocatable :: sfaces
        character(len=16), dimension(:), allocatable :: surfacenames

        allocate(surfacenames(meshdetails(7)))
        allocate(sfaces(meshdetails(2),4))
        sfaces = 0
        do i=1,meshdetails(7)
            read(unum,'(i8,1x,a)') surfnum, surfacenames(i)
            if(mod(surfnum,5) == 0) then
                countersize = surfnum/5
            else
                countersize = 1 + surfnum/5
            end if
            do j=1,countersize
                read(unum, '(5(2x,i8,1x,i1))') tempmat
                do k=1,5
                    if(tempmat(2*k-1) /= 0) then
                        sfaces(tempmat(2*k-1),tempmat(2*k)) = i
                    end if
                end do
            end do
        end do
    end subroutine readmeshsurfaces

	subroutine readSurfaces(unum,mDets,surfaces)
		integer :: i,j,k,numSurf,numFc,ct,cs,unum,mDets(7),temp(10)
		type(surfaceData),allocatable :: surfaces(:)

		numSurf = mDets(7)
		allocate(surfaces(mDets(7)))
		do i=1,numSurf
			surfaces%sfNumber = i
            read(unum,'(i8,1x,a)') numFc, surfaces(i)%sfName
			allocate(surfaces(i)%elNum(numFc))
			allocate(surfaces(i)%fcNum(numFc))
			if(mod(numFc,5) == 0) then
                cs = numFc/5
            else
                cs = 1 + numFc/5
            end if
			ct = 0
            do j=1,cs
                read(unum,'(5(2x,i8,1x,i1))') temp
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
