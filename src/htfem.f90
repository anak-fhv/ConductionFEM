module htfem

	use readmesh
	use element
	use boundary
	use assembly
	use solver
	use postproc
	use transient

	implicit none

	contains

	subroutine fem()
		integer,parameter :: resfilenum = 101, outfilenum = 102			! Files to store results and other output
		integer :: nNodes,nElems,nDoms,nSurfs,elDom,fcBytype,i,j,k,	&	! Prefixes: n=>number, el=>element, fc=>face
				   meshVals(7),elNodes(4),elByfaces(4),iter				! by=>boundary, sf=>surface, gn=>generation, no=>node
		integer,allocatable :: doElems(:),byCs(:),stRowPtr(:),		&	! do=>domain, sy=>global system
							   stCols(:),cpRowPtr(:),cpCols(:),		&
							   connTab(:,:),sfElems(:,:)
		real(8),parameter :: kDefault = 1.d0
		real(8) :: elVol,tAmbient,gnVal,elK,tc,qBHigh,qBLow,		&
				   byTemp(4),bySrc(4),gnSrc(4),elVerts(4,3),		&
				   elSpfns(4,4),elSt(4,4),bySt(4,4),elCp(4,4)
		real(8),allocatable :: domKs(:),sfVals(:),sySt(:),sySrc(:), &
							   syTvals(:),noVerts(:,:),reVals(:),	&
							   vF(:),domRhos(:),domCs(:),syCp(:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"results.out",	&
								  outfile = objdir//"outputs.out"
		character(16),allocatable :: sfFcname(:)
		logical,parameter :: gnDefault = .false.,trDefault = .false.
		logical ::	gnUser,trUser,useRK
		type(noderow) :: noElemPart(4),cpElemPart(4)
		type(noderow),allocatable :: stNo(:),cpNo(:)

		call getmeshdata(meshVals,noVerts,connTab,doElems,			&
		sfFcname,sfElems)

		write(*,'(a)') "Meshdetails received"

		nNodes = meshVals(1)
		nElems = meshVals(2)
		nDoms  = meshVals(6)
		nSurfs = meshVals(7)

		write(*,*) "System size: ", nNodes
		allocate(sySrc(nNodes))
		allocate(syTvals(nNodes))
		allocate(stNo(nNodes))

		if(nDoms .gt. 1) then
			allocate(vF(nDoms))
			vF = 0.d0
		end if

		sySrc = 0.d0
		syTvals = 0.d0

		call readboundaryconditions(meshVals,byCs,domKs,domRhos,	&
		domCs,sfVals,tAmbient,trUser,gnUser)

		trUser = .true.

		if(trUser) then
			allocate(cpNo(nNodes))
		end if

		do i=1,nElems

!	Get element data from mesh values
			elNodes = connTab(i,:)
			elDom = doElems(i)
			elVerts = noVerts(elNodes,:)
			elK = domKs(elDom)

!	Call stiffness functions
			call shapefunctions(elVerts,elVol,elSpfns)
			call elementstiffness(elSpfns,elVol,elK,elSt)

			if(trUser) then
				call getcapacitance(elDom,elVol,domCs,domRhos,elCp)
			end if

			if(nDoms .gt. 1) then
				vF(elDom) = vF(elDom) + elVol
			end if

!	Check boundary conditions
			elByfaces = sfElems(i,:)
			if(all(elByfaces==0)) then
				continue
			else
				call boundaryconditions(elVerts,elByfaces,byCs,		&
				sfVals,tAmbient,bySt,bySrc,byTemp)
				elSt = elSt + bySt
				call addtoglobaltemperature(syTvals,elNodes,byTemp)
				call addtoglobalforce(sySrc,elNodes,bySrc)
			end if

!	Check generation sources
			if(gnUser) then
				call uniformgeneration(gnVal,elVol,gnSrc)
				call addtoglobalforce(sySrc,elNodes,gnSrc)
			end if

!	Assemble the noderows
			noElemPart = stNo(elNodes)
			call assemble_noderows(noElemPart,elNodes,elSt)
			stNo(elNodes) = noElemPart

			cpElemPart = cpNo(elNodes)
			call assemblecapacitance(cpElemPart,elNodes,elCp)
			cpNo(elNodes) = cpElemPart

		end do

		if(trUser) then
!			call readinitialvalues(syTvals)
			syTvals = 100.d0
		else
			call setupfinalequations(stNo,sySrc,syTvals)
		end if

		call collapse_noderows(stNo,sySt,stCols,stRowPtr)

		if(trUser) then
			call collapse_noderows(cpNo,syCp,cpCols,cpRowPtr)
		end if

		deallocate(stNo)
		deallocate(cpNo)

		write(*,*) "Entered solution step"

		if(trUser) then
			useRK = .false.
			call transientsolve(sySt,stRowPtr,stCols,syCp,cpRowPtr,		&
			cpCols,sySrc,useRK,syTvals,noVerts)
		else
			call bicgstab(sySt,stRowPtr,stCols,sySrc,100000,revals,iter)
			write(*,'(a,i5,2x,a)') "This program took: ",iter,			&
			"iterations to converge."

			open(resfilenum,file=resfile)

			do i=1,nNodes
				write(resfilenum,'(3(f9.4,2x),f9.4)')noVerts(i,1:3), 	&
				revals(i)
			end do

			write(resfilenum,*)

			call getflowrates(noVerts,connTab,doElems,domKs,sfElems,	&
			reVals,(/1,2/),(/3,4/),3,qBLow,qBHigh)

			if(nDoms .eq. 2) then
				write(resfilenum,'(a,f9.4)') "Sample porosity:",		&
				vF(1)/(sum(vF))
			end if

			write(resfilenum,*) "Fluxes:"
			write(resfilenum,'(a,f9.4)') "Boundary low:", qBLow
			write(resfilenum,'(a,f9.4)') "Boundary high:", qBHigh

			close(resfilenum)
		end if

	end subroutine fem

    subroutine getmeshdata(meshdetails,vertices,connectivity,		&
	domainelements,surfacenames,surfacefaces)
        integer,parameter :: unitnumber = 103
		integer,dimension(7) :: meshdetails
        integer,dimension(:,:),allocatable :: connectivity, surfacefaces
        integer,dimension(:),allocatable :: domainelements
        character(len=16),dimension(:),allocatable :: surfacenames
        real(8),dimension(:,:),allocatable :: vertices
        real(8),dimension(:),allocatable :: boundaryvalues

        call openmeshfile(unitnumber, 'a.msh')
        call readmeshdetails(unitnumber,meshdetails)
        call readmeshvertices(unitnumber,meshdetails, vertices)
        call readmeshconnectivity(unitnumber,meshdetails, 			&
		connectivity)
        call readmeshdomains(unitnumber,meshdetails,domainelements)
        call readmeshsurfaces(unitnumber,meshdetails,surfacefaces,	&
		surfacenames)
        call closemeshfile(unitnumber)
    end subroutine getmeshdata

	subroutine addtoglobaltemperature(Tvals,elnodes,btemp)
		real(8),dimension(:),intent(inout) :: Tvals
		real(8),dimension(4) :: btemp
		integer,dimension(4) :: elnodes

		if(any(btemp.ne.0.0d0)) then
			Tvals(elnodes) = btemp
		end if
	end subroutine addtoglobaltemperature

	subroutine addtoglobalforce(gF,elnodes,elf)
		real(8),dimension(:),intent(inout) :: gF
		real(8),dimension(4) :: elf
		integer,dimension(4) :: elnodes
		integer :: i

		do i=1,4
			gF(elnodes(i)) = gF(elnodes(i)) + elf(i)
		end do
	end subroutine addtoglobalforce

	subroutine uniformgeneration(gval,elvol,gcontrib)
		real(8) :: gval,elvol
		real(8),dimension(4) :: gcontrib

		gcontrib = gval*(elvol/24.0d0)*(/1,1,1,1/)
	end subroutine uniformgeneration

	subroutine setupfinalequations(stNo,sySrc,syTvals)
		integer :: i
		real(8) :: sySrc(:),syTvals(:)
		type(noderow):: stNo(:)

		do i=1,size(stNo,1)
			if(syTvals(i) .ne. 0.0d0) then
				if(allocated(stNo(i)%col)) then
					deallocate(stNo(i)%val)
					deallocate(stNo(i)%col)
				end if
				allocate(stNo(i)%col(1))
				allocate(stNo(i)%val(1))
				stNo(i)%col(1) = i
				stNo(i)%val(1) = 1.0d0
				sySrc(i) = syTvals(i)
			end if
		end do
	end subroutine setupfinalequations

end module htfem
