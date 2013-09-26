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
        real :: boundaryval, facearea, elementvolume, tambient, generationrate
        real,dimension(:,:),allocatable :: vertices, elementcoords, elementgradient,  &
        conductivestiffness, convectivestiffness, elementstiffness, globalstiffness, stiffinv
        real,dimension(:,:,:),allocatable :: conductivities,stels
        real,dimension(:),allocatable :: boundaryvalues, elementforce, &
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
        allocate(boundaryconditions(numbd))
        allocate(boundaryvalues(numbd))
        allocate(temperature(numno))
        allocate(temp2(numno))
        allocate(globalflux(numno))
        allocate(conductivities(numdo,3,3))
        allocate(stels(numel,4,4))

!        globalstiffness = 0.0d0
        globalflux = 0.0d0
        temperature = 0.0d0

        call readboundaryconditions(meshdetails, boundaryconditions, &
        conductivities, boundaryvalues, tambient, generation, generationrate)
        print *, 'Boundaryconditions: ', boundaryconditions
        print *, 'Surfacefaces: ', surfacefaces(207,:)
        do i=1,3
            print *, 'Element: ', i
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
            !call assembleglobalstiffness(elementstiffness, elementnodes, globalstiffness)
            do co=1,4
                do bo=1,4
                    m = elementnodes(co)
                    n = elementnodes(bo)
!                    print *, 'm,n: ', m,n
                    globalstiffness(m,n) = globalstiffness(m,n) + elementstiffness(co,bo)
                end do
            end do
            if(generation) then
                call getelementflux(elementgradient,generationrate,elementvolume,elementflux)
                call assembleelementfluxes(globalflux,elementflux,elementnodes)
            end if
        end do
        co = 0
        bo = 0
        do i= 1,numno
            if(globalflux(i) == 0.0d0) then
                co = co+1
            end if
            if(abs(vertices(i,3)) > 0.9999d0) then
                bo = bo+1
            end if
        end do
    end subroutine finitesolution

    subroutine getmeshdata(meshdetails,vertices,connectivity,domainelements,surfacenames,surfacefaces)
        integer, parameter :: unitnumber = 111
        integer, dimension(:,:), allocatable :: connectivity, surfacefaces
        integer, dimension(:), allocatable :: domainelements, meshdetails
        character(len=16), dimension(:), allocatable :: surfacenames
        real, dimension(:,:), allocatable :: vertices
        real, dimension(:), allocatable :: boundaryvalues

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
        character(len=*),parameter :: filename='/home/anak/Documents/Old_Data/omg/trials/datafile.dat'
        integer,dimension(7) :: meshdetails
        integer :: numdomains,numboundaries,i,j
        logical :: generation
        integer,dimension(:) :: boundaryconditions
        real :: tambient, generationrate
        real,dimension(:) :: boundaryvalues
        real,dimension(:,:,:) :: conductivities

        numdomains = meshdetails(6)
        numboundaries = meshdetails(7)
!        allocate(conductivities(numdomains,3,3))
!        allocate(boundaryconditions(numboundaries))
!        allocate(boundaryvalues(numboundaries))
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

    subroutine getboundaryconditions(surfacenames, boundaryconditions, boundaryvalues)
        character(len=16), dimension(:) :: surfacenames
        integer :: i
        integer, dimension(:), allocatable :: boundaryconditions
        real, dimension(:), allocatable :: boundaryvalues

        allocate(boundaryconditions(size(surfacenames,1)))
        allocate(boundaryvalues(size(surfacenames,1)))
        print *, 'Please enter the type of boundary that each surface is.'
        print *, 'Use 1 for a prescribed temperature boundary, 2 for a ', &
        'prescribed flux boundary and 3 for a convective transfer boundary'
        print *, 'Use 4 if this is an interface within the material.'
        do i=1,size(surfacenames,1)
            boundaryconditions(i) = getuserboundarycondition(surfacenames(i))
            if(boundaryconditions(i) == 4) then
                boundaryvalues(i) = 0.0
            else
                boundaryvalues(i) = getuserboundaryvalue(surfacenames(i),boundaryconditions(i))
            end if
        end do
    end subroutine getboundaryconditions

!This is a recursive function in order to force the user to choose the right
!values for boundary condition types. Here, we parametrically use 1 for a
!prescribed temperature BC, 2 for a prescribed flux and 3 for a prescribed
!heat transfer coefficient. The BC 4 refers to an interface between phases.

    recursive function getuserboundarycondition(surface) result(condition)
        character(len=16) :: surface
        integer :: condition
        integer, dimension(4), parameter :: allowed_values = (/1,2,3,4/)

        print *, 'Enter the type of boundary that ', trim(surface), ' is: '
        read(*,'(i1)') condition
        print *, 'Condition: ', condition
        if(condition<1 .or. condition>4) then
            print *, 'Error, type specified unrecognized. Retrying...'
            condition = getuserboundarycondition(surface)
        end if
        return
    end function getuserboundarycondition

!This is a recursive function because we need to force the user to enter +ve
!values of temperature.

    recursive function getuserboundaryvalue(surface,condition) result(value)
        character(len=16) :: surface
        integer :: condition
        real :: value

        if(condition == 1) then
            print *, 'Enter temperature at the boundary ', trim(surface)
            read(*,*) value
            if(value < 0.0) then
                print *, 'Cannot prescribe a negative temperature. Retrying...'
                value = getuserboundaryvalue(surface, condition)
            end if
        elseif(condition == 2) then
            print *, 'Enter the value of flux at the boundary ', trim(surface), ' including the sign'
            read(*,*) value
        elseif(condition == 3) then
            print *, 'Enter the value of the heat transfer coefficient at the boundary ', trim(surface)
            read(*,*) value
        end if
        return
    end function getuserboundaryvalue

    subroutine getconductivities(numdomains,conductivities)
        integer :: numdomains, i,j
        real, dimension(:,:,:), allocatable :: conductivities

        allocate(conductivities(numdomains,3,3))
        do i=1,numdomains
            do j=1,3
                print *, 'Enter the conductivity row ', j, ' of the domain ', i
                read(*,*) conductivities(i,j,:)
            end do
        end do

    end subroutine getconductivities



    subroutine assembletemperature(temperature,facenumber,elementnodes,boundaryvalue)
        real,dimension(:),intent(inout) :: temperature
        real :: boundaryvalue
        integer,dimension(4) :: elementnodes
        integer,dimension(3) :: facenodes
        integer :: facenumber

        facenodes = getfacenodeslocal(facenumber)
        temperature(elementnodes(facenodes)) = boundaryvalue
    end subroutine assembletemperature

    subroutine getelementflux(elementgradient,boundaryval,facearea,elementflux)
        real,dimension(4,4) :: elementgradient
        real,dimension(4) :: elementflux
        real :: boundaryval,facearea
!        print *, 'Elementgradient: ', elementgradient(1,:)
!        print *, 'boundaryval: ', boundaryval
!        print *, 'Facearea: ', facearea
        elementflux = boundaryval*facearea*elementgradient(1,:)
!        print *, 'Elementflux: ', elementflux
    end subroutine getelementflux

    subroutine assembleelementfluxes(globalflux,elementflux,elementnodes)
        real,dimension(:),intent(inout) :: globalflux
        real,dimension(4) :: elementflux
        integer,dimension(4) :: elementnodes
        integer :: i
!        print *, 'GF: ', size(globalflux,1)
        do i=1,4
!            print *, globalflux(elementnodes(i))
            globalflux(elementnodes(i)) = globalflux(elementnodes(i)) + elementflux(i)
        end do
    end subroutine assembleelementfluxes

    function getboundaryarea(coords,facenumber,faceloc) result(area)
        real,dimension(4,3) :: coords
        real,dimension(3,3):: facecoords
        integer :: boundary, faceloc, facenumber
        integer,dimension(3) :: facenodes
        real :: area
        
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
        real, dimension(4,3) :: elementcoords
        real, dimension(4,4) :: volumematrix, elementgradient
        real, dimension(3,3) :: minor
        real :: elementvolume
        integer :: i,j
        integer,dimension(2,1) :: trial

        trial = reshape((/2,3/),(/2,1/))
        elementvolume = 0.0d0
        volumematrix = 0.0d0
        volumematrix(1:4,2:4) = elementcoords
        volumematrix(:,1) = 1.0d0
        print *, 'volumematrix: '
        print *, volumematrix(1,:)
        print *, volumematrix(2,:)
        print *, volumematrix(3,:)
        print *, volumematrix(4,:)
        call solvedeterminant(volumematrix, elementvolume)
        elementvolume = abs(elementvolume/6.0d0)
!        call gaussjordaninverse(volumematrix,elementgradient)
        call getelementgradient(elementcoords,elementgradient)
        elementgradient = transpose(elementgradient)/(6.0d0*elementvolume)
        print *, 'Elementgradient: '
        print *, elementgradient(1,:)
        print *, elementgradient(2,:)
        print *, elementgradient(3,:)
        print *, elementgradient(4,:)
    end subroutine getelementbasisdata

    subroutine getconductivestiffness(elementgradient,conductivity, conductivestiffness)
        real, dimension(4,4) :: elementgradient, conductivestiffness
        real, dimension(4,3) :: intermediate
        real, dimension(3,3) :: conductivity

        intermediate = matmul(transpose(elementgradient(2:4,:)),conductivity)
        conductivestiffness = matmul(intermediate,elementgradient(2:4,:))
    end subroutine getconductivestiffness

    subroutine getconvectivestiffness(elementgradient,area,hvalue,convectivestiffness)
        real, dimension(4,3) :: elementcoords
        real, dimension(4,4) :: elementgradient, convectivestiffness
        real :: hvalue, area
        integer :: j

        convectivestiffness = hvalue*area*matmul(transpose(elementgradient(1:1,:)),elementgradient(1:1,:))
    end subroutine getconvectivestiffness

    subroutine assembleglobalstiffness(elementstiffness, elementnodes, globalstiffness)
        real, dimension(4,4) :: elementstiffness
        real, dimension(:,:), intent(inout) :: globalstiffness
        integer, dimension(4) :: elementnodes
        integer :: i,j,m,n

        print *, 'Globalstiffness: '
        print *, globalstiffness(elementnodes(1),elementnodes(1))
        print *, globalstiffness(elementnodes(4),elementnodes(4))
        do i=1,4
            m = elementnodes(i)
            do j=1,4
                n = elementnodes(j)
                globalstiffness(m,n) = globalstiffness(m,n) + elementstiffness(i,j)
            end do
        end do
        print *, 'Globalstiffness modified?'
        print *, globalstiffness(elementnodes(1),elementnodes(1))
        print *, globalstiffness(elementnodes(4),elementnodes(4))
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

    subroutine getelementgradient(E,G)
        real,dimension(4,3) :: E
        real,dimension(4,4) :: G

        G(1,1) = E(2,1)*(E(3,2)*E(4,3)-E(4,2)*E(3,3)) - E(2,2)*(E(3,1)*E(4,3)-E(4,1)*E(3,3)) + E(2,3)*(E(3,1)*E(4,2)-E(4,1)*E(3,2))
        G(1,2) = E(3,1)*(E(4,2)*E(1,3)-E(1,2)*E(4,3)) - E(3,2)*(E(4,1)*E(1,3)-E(1,1)*E(4,3)) + E(3,3)*(E(4,1)*E(1,2)-E(1,1)*E(4,2))
        G(1,3) = E(4,1)*(E(1,2)*E(2,3)-E(2,2)*E(1,3)) - E(4,2)*(E(1,1)*E(2,3)-E(2,1)*E(1,3)) + E(4,3)*(E(1,1)*E(2,2)-E(2,1)*E(1,2))
        G(1,4) = E(1,1)*(E(2,2)*E(3,3)-E(3,2)*E(2,3)) - E(1,2)*(E(2,1)*E(3,3)-E(3,1)*E(2,3)) + E(1,3)*(E(2,1)*E(3,2)-E(3,1)*E(2,2))

        G(2,1) = (E(2,2)-E(4,2))*(E(3,3)-E(4,3)) - (E(3,2)-E(4,2))*(E(2,3)-E(4,3))
        G(2,2) = (E(3,2)-E(4,2))*(E(1,3)-E(4,3)) - (E(1,2)-E(4,2))*(E(3,3)-E(4,3))
        G(2,3) = (E(1,2)-E(4,2))*(E(2,3)-E(4,3)) - (E(2,2)-E(4,2))*(E(1,3)-E(4,3))
        G(2,4) = (G(2,1) + G(2,2) + G(2,3))

        G(3,1) = (E(3,1)-E(4,1))*(E(2,3)-E(4,3)) - (E(2,1)-E(4,1))*(E(3,3)-E(4,3))
        G(3,2) = (E(1,1)-E(4,1))*(E(3,3)-E(4,3)) - (E(3,1)-E(4,1))*(E(1,3)-E(4,3))
        G(3,3) = (E(2,1)-E(4,1))*(E(1,3)-E(4,3)) - (E(1,1)-E(4,1))*(E(2,3)-E(4,3))
        G(3,4) = -(G(3,1) + G(3,2) + G(3,3))

        G(4,1) = (E(2,1)-E(4,1))*(E(3,2)-E(4,2)) - (E(3,1)-E(4,1))*(E(2,2)-E(4,2))
        G(4,2) = (E(3,1)-E(4,1))*(E(1,2)-E(4,2)) - (E(1,1)-E(4,1))*(E(3,2)-E(4,2))
        G(4,3) = (E(1,1)-E(4,1))*(E(2,2)-E(4,2)) - (E(2,1)-E(4,1))*(E(1,2)-E(4,2))
        G(4,4) = -(G(4,1) + G(4,2) + G(4,3))
    end subroutine getelementgradient

    subroutine getelementforcecomponent()
    end subroutine getelementforcecomponent

end module finitesolver
