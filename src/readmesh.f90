module readmesh
    use numericalmethods
    implicit none
    contains

    subroutine openmeshfile(unitnumber, filename)
        integer :: unitnumber
        character(len = *) :: filename

        open(unit=unitnumber, file=filename, form='formatted', status='old')
    end subroutine openmeshfile

    subroutine closemeshfile(unitnumber)
        integer :: unitnumber

        close(unitnumber)
    end subroutine closemeshfile

    subroutine readmeshdetails(unitnumber, meshdetails)
        integer :: unitnumber
        integer, dimension(7) :: meshdetails

        read(unitnumber, *)
        read(unitnumber, *)
        read(unitnumber, *) meshdetails
    end subroutine readmeshdetails

    subroutine readmeshvertices(unitnumber, meshdetails, vertices)
        integer :: unitnumber
        integer, dimension(7) :: meshdetails
        integer :: i
        double precision, dimension(:,:), allocatable :: vertices

        allocate(vertices(meshdetails(1), 3))
        do i=1,meshdetails(1)
            read(unitnumber, *) vertices(i,:)
        end do
    end subroutine readmeshvertices

! This subroutine has been written only for a tetrahedral element mesh

    subroutine readmeshconnectivity(unitnumber, meshdetails, connectivity)
        integer :: unitnumber,i
        integer, dimension(7) :: meshdetails
        integer, dimension(:,:), allocatable :: connectivity

        allocate(connectivity(meshdetails(2), 4))
        do i=1,meshdetails(2)
            read(unitnumber, *) connectivity(i,:)
        end do
    end subroutine readmeshconnectivity

! This subroutine simply notes on each element of the mesh the domain to
! which it belongs. So it is an n x 1 array.

    subroutine readmeshdomains(unitnumber, meshdetails, domainelements)
        integer :: unitnumber,domainsize,countersize,i,j,k
        integer, dimension(7) :: meshdetails
        integer, dimension(10) :: tempmat
        integer, dimension(:), allocatable :: domainelements

        allocate(domainelements(meshdetails(2)))
        tempmat = 0
        do i=1,meshdetails(6)
            read(unitnumber,'(i8)') domainsize
            if(mod(domainsize,10) == 0) then
                countersize = domainsize/10
            else
                countersize = 1 + (domainsize/10)
            end if
            do j=1,countersize
                read(unitnumber,'(10(1x,i8))') tempmat
                do k=1,10
                    if(tempmat(k) /= 0) then
                        domainelements(tempmat(k)) = i
                    end if
                end do
            end do
        end do
    end subroutine readmeshdomains

! This subroutine has been written only for tetrahedral elements. It uses
! an n x 4 array, and on each ith member of this array, for each of the 4
! tetrahedron faces, it notes a 0 if the face is not on a boundary or the
! boundary number (as sequentially occurs) if it is on a boundary

    subroutine readmeshsurfaces(unitnumber, meshdetails, surfacefaces, surfacenames)
        integer :: unitnumber,surfnum,countersize,i,j,k, currentsize, newsize, oldsize
        integer, dimension(7) :: meshdetails
        integer, dimension(10) :: tempmat
        integer, dimension(:,:), allocatable :: surfacefaces
        character(len=16), dimension(:), allocatable :: surfacenames

        allocate(surfacenames(meshdetails(7)))
        allocate(surfacefaces(meshdetails(2),4))
        surfacefaces = 0
        do i=1,meshdetails(7)
            read(unitnumber,'(i8,1x,a)') surfnum, surfacenames(i)
            if(mod(surfnum,5) == 0) then
                countersize = surfnum/5
            else
                countersize = 1 + surfnum/5
            end if
            do j=1,countersize
                read(unitnumber, '(5(2x,i8,1x,i1))') tempmat
                do k=1,5
                    if(tempmat(2*k-1) /= 0) then
                        surfacefaces(tempmat(2*k-1),tempmat(2*k)) = i
                    end if
                end do
            end do
        end do
    end subroutine readmeshsurfaces

end module readmesh
