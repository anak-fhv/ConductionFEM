module boundary

	use element

	implicit none
	contains

!-------------------------------------------------------------------
! Routine to read boundary conditions from the data file
!-------------------------------------------------------------------
    subroutine readboundaryconditions(meshdetails, boundaryconditions, &
    kvals, boundaryvalues, Tambient, generation, generationrate)
        integer,parameter :: datafilenum=222
        character(len=*),parameter :: filename='datafile.dat'
        integer,dimension(7) :: meshdetails
        integer :: numdomains,numboundaries,i,j
        logical :: generation
        integer,dimension(:),allocatable :: boundaryconditions
        real(8) :: Tambient
		real(8),optional :: generationrate
        real(8),dimension(:),allocatable :: kvals,boundaryvalues

        numdomains = meshdetails(6)
        numboundaries = meshdetails(7)
		allocate(boundaryconditions(numboundaries))
		allocate(kvals(numdomains))
		allocate(boundaryvalues(numboundaries))
        boundaryvalues = 0.0d0
        boundaryconditions = 4
        open(datafilenum,file=filename,status='old')
        read(datafilenum,*)
        read(datafilenum,*)
        do i=1,numdomains
			read(datafilenum,*)
            read(datafilenum,*) kvals(i)
        end do
		read(datafilenum,*)
        read(datafilenum,*) boundaryconditions
        read(datafilenum,*)
        read(datafilenum,*) boundaryvalues
        read(datafilenum,*)
        read(datafilenum,*) Tambient
        read(datafilenum,*)
        read(datafilenum,*) generation
        if(generation) then
            read(datafilenum,*)
            read(datafilenum,*) generationrate
        end if
    end subroutine readboundaryconditions
!-------------------------------------------------------------------
! END subroutine to read boundary conditions from data file
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! Giant subroutine for the handling of boundary conditions
!-------------------------------------------------------------------
	subroutine boundaryconditions(ec,bofcs,bcs,bvals,Tamb,bstiff,		&
	bforce,btemp)
		real(8),dimension(4,3) :: ec
		real(8),dimension(4,4) :: bstiff,elht
		real(8),dimension(4) :: bforce,btemp,elq,hta
		real(8),dimension(:),intent(in) :: bvals
		real(8) :: Tamb,bv
		integer,dimension(:),intent(in) :: bcs
		integer,dimension(4) :: bofcs
		integer,dimension(3) :: fcnodes
		integer :: i,j,bloc,fbtype,n1,n2

		bstiff = 0.0d0
		bforce = 0.0d0
		btemp = 0.0d0
		elq = 0.0d0
		hta = 0.0d0
		elht = 0.0d0
		do i=1,4
			if(bofcs(i) /= 0) then
				bloc = bofcs(i)
				fbtype = bcs(bloc)
				bv = bvals(bloc)
				call bfacenodes(i,fcnodes)
				if(bv.ne.0.0d0) then
					if(fbtype == 1) then
						btemp(fcnodes) = bv
					elseif(fbtype == 2) then
						call fluxboundary(ec,fcnodes,bv,elq)
						bforce = bforce - elq
					elseif(fbtype == 3) then
						call convectiveboundary(ec,fcnodes,bv,Tamb,hta,elht)
						bforce = bforce + hta
						bstiff = bstiff + elht
					else
						continue
					end if
				end if
			else
				continue
			end if
		end do
	end subroutine boundaryconditions
!-------------------------------------------------------------------
! END subroutine for the handling of boundary conditions
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! Smaller subroutines called from within the giant routine
!-------------------------------------------------------------------
	subroutine temperatureboundary(fcnum,Telem)
		real(8),dimension(4) :: Telem
		integer :: fcnum
	end subroutine temperatureboundary

	subroutine fluxboundary(ec,fcnodes,bv,elq)
		real(8),dimension(4,3) :: ec
		real(8),dimension(3,3) :: fc
		real(8),dimension(4) :: elq
		real(8) :: fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: i,j,k,el

		elq = 0.0d0
		fc = ec(fcnodes,:)
		fcarea = facearea(fc)
		elq(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv
	end subroutine fluxboundary

	subroutine convectiveboundary(ec,fcnodes,bv,Tamb,hta,elht)
		real(8),dimension(4,3) :: ec
		real(8),dimension(4,4) :: sp,elht,surfint
		real(8),dimension(3,3) :: fc
		real(8),dimension(4) :: hta
		real(8) :: Tamb,fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: fnum,i,bl

		hta = 0.0d0
		elht= 0.0d0
		fc = ec(fcnodes,:)
		do i=1,3
		end do
		fcarea = facearea(fc)
		hta(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv*Tamb
		surfint = shapefuncsquaresurfint(fcnodes)
		elht = (2.0d0*fcarea/24.0d0)*bv*surfint
	end subroutine convectiveboundary
!-------------------------------------------------------------------
! END smaller subroutines
!-------------------------------------------------------------------

end module boundary
