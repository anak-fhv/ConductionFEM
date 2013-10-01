module finitesolver
    use readmesh
    use numericalmethods
    implicit none
    contains

    subroutine finitesolution()
        integer,dimension(:,:),allocatable :: connectivity, surfacefaces
        integer,dimension(:),allocatable :: meshdetails,domainelements,elementnodes, &
        boundaryconditions, facenodes
        integer :: i,j,k, numno, numel, numdo, numbd, boundary, boundarytype, faceloc, co, bo, m,n
        double precision :: boundaryval, facearea, elementvolume, tambient, generationrate
        double precision,dimension(:,:),allocatable :: vertices, elementcoords, elementgradient,  &
        conductivestiffness, convectivestiffness, elementstiffness, globalstiffness, stiffinv
        double precision,dimension(:,:,:),allocatable :: conductivities
        double precision,dimension(:),allocatable :: boundaryvalues, elementforce, &
        elementflux,globalflux,temperature,temp2
        character(len=16), dimension(:), allocatable :: surfacenames
        logical :: generation

        print *, 'Entered finite solver'
        call getmeshdata(meshdetails,vertices,connectivity,domainelements,surfacenames,surfacefaces)
        print *, 'Meshdetails: ', meshdetails
        numno = meshdetails(1)
        numel = meshdetails(2)
        numdo = meshdetails(6)
        numbd = meshdetails(7)
        allocate(elementcoords(4,3))
        allocate(elementgradient(4,4))
        allocate(elementstiffness(4,4))
        allocate(conductivestiffness(4,4))
        allocate(convectivestiffness(4,4))
        allocate(globalstiffness(numno,numno))
        allocate(stiffinv(numno,numno))
        allocate(elementnodes(4))
        allocate(elementflux(4))
        allocate(temperature(numno))
        allocate(temp2(numno))
        allocate(globalflux(numno))

        globalstiffness = 0.0d0
        globalflux = 0.0d0
        temperature = 0.0d0

        call readboundaryconditions(meshdetails, boundaryconditions, &
        conductivities, boundaryvalues, tambient, generation, generationrate)

        print *, 'Boundaryconditions: ', boundaryconditions

        do i=1,numel
            elementstiffness = 0.0d0
            conductivestiffness = 0.0d0
            convectivestiffness = 0.0d0
            elementnodes = connectivity(i,:)
            elementcoords = vertices(elementnodes,:)
            call getelementbasisdata(elementcoords,elementgradient,elementvolume)
            call getconductivestiffness(elementgradient,conductivities,conductivestiffness)
            do j=1,4
                boundary = surfacefaces(i,j)
                if(boundary /= 0) then
                    boundarytype = boundaryconditions(boundary)
                    boundaryval = boundaryvalues(boundary)
                    if(boundarytype == 1) then
                        call assembletemperature(temperature,j,elementnodes,boundaryval)
                    end if
                    if(boundarytype == 2 .or. boundarytype == 3) then
                        faceloc = getfacelocationindex(surfacefaces(i,j),surfacenames)
                        facearea = getboundaryarea(elementcoords,j,faceloc)
                        if(boundarytype == 2) then
                            call getelementflux(elementgradient,boundaryval,facearea,elementflux)
                            call assembleelementfluxes(globalflux,elementflux,elementnodes)
                        else
                            call getconvectivestiffness(elementgradient,boundaryval,facearea,convectivestiffness)
                            call getelementflux(elementgradient,boundaryval,facearea,elementflux)
                            elementflux = tambient*elementflux
                            call assembleelementfluxes(globalflux,elementflux,elementnodes)
                        end if
                    end if
                    if(boundarytype == 4) then
                        continue
                    end if
                end if
            end do
            elementstiffness = conductivestiffness + convectivestiffness
            call assembleglobalstiffness(elementstiffness, elementnodes, globalstiffness)
            if(generation) then
                call getelementflux(elementgradient,generationrate,elementvolume,elementflux)
                call assembleelementfluxes(globalflux,elementflux,elementnodes)
            end if
        end do
        call solvefinalequations(globalstiffness,globalflux,temperature)
    end subroutine finitesolution

    subroutine solvefinalequations(globalstiffness,globalflux,temperature)
        double precision, dimension(:,:) :: globalstiffness
        double precision, dimension(:) :: globalflux, temperature
        integer :: m, n, info, lda, ldb, nrhs
        integer, dimension(:), allocatable :: ipiv
        character :: trans='N'

        m = size(globalstiffness,1)
        n = size(globalstiffness,2)
        lda = m
        ldb = m
        nrhs = 1
        allocate(ipiv(m))
        call dgesv(m, nrhs, globalstiffness, lda, ipiv, globalflux, ldb, info)
        print *, 'Info: ', info
        temperature = globalflux
    end subroutine solvefinalequations

    subroutine getmeshdata(meshdetails,vertices,connectivity,domainelements,surfacenames,surfacefaces)
        integer, parameter :: unitnumber = 111
        integer, dimension(:,:), allocatable :: connectivity, surfacefaces
        integer, dimension(:), allocatable :: domainelements, meshdetails
        character(len=16), dimension(:), allocatable :: surfacenames
        double precision, dimension(:,:), allocatable :: vertices
        double precision, dimension(:), allocatable :: boundaryvalues

        allocate(meshdetails(7))
        call openmeshfile(unitnumber, '/home/anak/Documents/Old_Data/omg/trials/a.msh')
        call readmeshdetails(unitnumber, meshdetails)
        call readmeshvertices(unitnumber, meshdetails, vertices)
        call readmeshconnectivity(unitnumber, meshdetails, connectivity)
        call readmeshdomains(unitnumber, meshdetails, domainelements)
        call readmeshsurfaces(unitnumber,meshdetails,surfacefaces,surfacenames)
        call closemeshfile(unitnumber)
    end subroutine getmeshdata

    subroutine readboundaryconditions(meshdetails, boundaryconditions, &
    conductivities, boundaryvalues, tambient, generation, generationrate)
        integer,parameter :: datafilenum=222
        character(len=*),parameter :: filename='datafile.dat'
        integer,dimension(7) :: meshdetails
        integer :: numdomains,numboundaries,i,j
        logical :: generation
        integer,dimension(:),allocatable :: boundaryconditions
        double precision :: tambient, generationrate
        double precision,dimension(:),allocatable :: boundaryvalues
        double precision,dimension(:,:,:),allocatable :: conductivities

        numdomains = meshdetails(6)
        numboundaries = meshdetails(7)
        allocate(conductivities(numdomains,3,3))
        allocate(boundaryconditions(numboundaries))
        allocate(boundaryvalues(numboundaries))
        conductivities = 0
        boundaryvalues = 0.0d0
        boundaryconditions = 4
        open(datafilenum,file=filename,status='old')
        read(datafilenum,*)
        read(datafilenum,*)
        do i=1,numdomains
            read(datafilenum,*)
            do j=1,3
                read(datafilenum,*) conductivities(i,j,:)
            end do
            read(datafilenum,*)
        end do
        read(datafilenum,*)
        read(datafilenum,*) boundaryconditions
        read(datafilenum,*)
        read(datafilenum,*) boundaryvalues
        read(datafilenum,*)
        read(datafilenum,*) tambient
        read(datafilenum,*)
        read(datafilenum,*) generation
        if(generation) then
            read(datafilenum,*)
            read(datafilenum,*) generationrate
        end if
    end subroutine readboundaryconditions

    subroutine assembletemperature(temperature,facenumber,elementnodes,boundaryvalue)
        double precision,dimension(:),intent(inout) :: temperature
        double precision :: boundaryvalue
        integer,dimension(4) :: elementnodes
        integer,dimension(3) :: facenodes
        integer :: facenumber

        facenodes = getfacenodeslocal(facenumber)
        temperature(elementnodes(facenodes)) = boundaryvalue
    end subroutine assembletemperature

    subroutine getelementflux(elementgradient,boundaryval,facearea,elementflux)
        double precision,dimension(4,4) :: elementgradient
        double precision,dimension(4) :: elementflux
        double precision :: boundaryval,facearea
        elementflux = boundaryval*facearea*elementgradient(1,:)
    end subroutine getelementflux

    subroutine assembleelementfluxes(globalflux,elementflux,elementnodes)
        double precision,dimension(:),intent(inout) :: globalflux
        double precision,dimension(4) :: elementflux
        integer,dimension(4) :: elementnodes
        integer :: i
        do i=1,4
            globalflux(elementnodes(i)) = globalflux(elementnodes(i)) + elementflux(i)
        end do
    end subroutine assembleelementfluxes

    function getboundaryarea(coords,facenumber,faceloc) result(area)
        double precision,dimension(4,3) :: coords
        double precision,dimension(3,3):: facecoords
        integer :: boundary, faceloc, facenumber
        integer,dimension(3) :: facenodes
        double precision :: area
        
        if(faceloc == 4) then
            print *, 'Face not recognized'
            area = 0
            return
        end if
        facenodes = getfacenodeslocal(facenumber)
        facecoords = coords(facenodes,:)
        facecoords(:,faceloc) = 1
        call solvedeterminant(facecoords,area)
        area = abs(area)
    end function getboundaryarea

    function getfacenodeslocal(facenumber) result(facenodes)
	integer :: facenumber
	integer,dimension(3) :: facenodes

        if(facenumber == 1) then
            facenodes = (/1,2,3/)
        elseif(facenumber == 2) then
            facenodes = (/1,2,4/)
        elseif(facenumber == 3) then
            facenodes = (/2,3,4/)
        else
            facenodes = (/4,1,3/)
        end if
    end function getfacenodeslocal

    subroutine getelementbasisdata(elementcoords, elementgradient, elementvolume)
        integer,parameter :: m=4, n=4, lda=4, lwork=256
        integer :: info
        integer, dimension(m) :: ipiv
        double precision, dimension(m) :: work
        double precision, dimension(4,3) :: elementcoords
        double precision, dimension(4,4) :: volumematrix, elementgradient
        double precision, dimension(3,3) :: minor
        double precision :: elementvolume

        volumematrix(2:4,1:4) = transpose(elementcoords)
        volumematrix(1,:) = 1.0d0
        call dgetrf(m,n,volumematrix,lda,ipiv,info)
        call dgetri(n,volumematrix,lda,ipiv,work,lwork,info)
        elementgradient = volumematrix
    end subroutine getelementbasisdata

    subroutine getconductivestiffness(elementgradient,conductivity, conductivestiffness)
        double precision, dimension(4,4) :: elementgradient, conductivestiffness
        double precision, dimension(4,3) :: intermediate
        double precision, dimension(3,3) :: conductivity

        intermediate = matmul(transpose(elementgradient(2:4,:)),conductivity)
        conductivestiffness = matmul(intermediate,elementgradient(2:4,:))
    end subroutine getconductivestiffness

    subroutine getconvectivestiffness(elementgradient,area,hvalue,convectivestiffness)
        double precision, dimension(4,3) :: elementcoords
        double precision, dimension(4,4) :: elementgradient, convectivestiffness
        double precision :: hvalue, area
        integer :: j

        convectivestiffness = hvalue*area*matmul(transpose(elementgradient(1:1,:)),elementgradient(1:1,:))
    end subroutine getconvectivestiffness

    subroutine assembleglobalstiffness(elementstiffness, elementnodes, globalstiffness)
        double precision, dimension(4,4) :: elementstiffness
        double precision, dimension(:,:), intent(inout) :: globalstiffness
        integer, dimension(4) :: elementnodes
        integer :: i,j,m,n

        do i=1,4
            m = elementnodes(i)
            do j=1,4
                n = elementnodes(j)
                globalstiffness(m,n) = globalstiffness(m,n) + elementstiffness(i,j)
            end do
        end do
    end subroutine assembleglobalstiffness

    function getfacelocationindex(bnum,surfacenames) result(faceloc)
        integer :: bnum, faceloc
        character(len=16),dimension(:) :: surfacenames

        if(index(surfacenames(bnum),'x') == 1) then
            faceloc = 1
        elseif(index(surfacenames(bnum),'y') == 1) then
            faceloc = 2
        elseif(index(surfacenames(bnum),'z') == 1) then
            faceloc = 3
        else
            faceloc = 4
        end if
    end function getfacelocationindex


end module finitesolver
