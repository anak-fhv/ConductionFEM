module htfem

	use readmesh
	use element
	use boundary
	use assembly
	use solver

	implicit none

	contains

	subroutine fem()
		integer,parameter :: resfilenum = 111, outfilenum = 222			! Files to store results and other output
		integer :: nNodes,nElems,nDoms,nSurfs,elDom,fcBytype,i,j,k,	&	! Prefixes: n=>number, el=>element, fc=>face
				   meshVals(7),elNodes(4),elByfaces(4),iter				! by=>boundary, sf=>surface, gn=>generation, no=>node
		integer,allocatable :: doElems(:),byCs(:),syRowPtr(:),		&	! do=>domain, sy=>global system
							   syCols(:),connTab(:,:),sfElems(:,:)
		real(8),parameter :: kDefault = 1.d0
		real(8) :: elVol,tAmbient,gnVal,elK,tc,byTemp(4),bySrc(4),	&
				   gnSrc(4),elVerts(4,3),elSpfns(4,4),elSt(4,4),	&
				   bySt(4,4)
		real(8),allocatable :: domKs(:),sfVals(:),sySt(:),sySrc(:), &
							   syTvals(:),noVerts(:,:),reVals(:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"results.out",	&
								  outfile = objdir//"outputs.out"
		character(16),allocatable :: sfFcname(:)
		logical,parameter :: gnDefault = .false.
		logical ::	gnUser
		type(noderow) :: noElemPart(4)
		type(noderow),allocatable :: noSys(:)

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
		allocate(noSys(nNodes))

		sySrc = 0.d0
		syTvals = 0.d0

		call readboundaryconditions(meshVals,byCs,domKs,sfVals,		&
		tAmbient,gnUser)

		do i=1,nElems

!	Get element data from mesh values
			elNodes = connTab(i,:)
			elDom = doElems(i)
			elVerts = noVerts(elNodes,:)
			elK = domKs(elDom)

!	Call stiffness functions
			call shapefunctions(elVerts,elVol,elSpfns)
			call elementstiffness(elSpfns,elVol,elK,elSt)

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
			noElemPart = noSys(elNodes)
			call assemble_noderows(noElemPart,elNodes,elSt)
			noSys(elNodes) = noElemPart

		end do

		deallocate(connTab)
		deallocate(doElems)

		call setupfinalequations(noSys,sySrc,syTvals)

		call collapse_noderows(noSys,sySt,syCols,syRowPtr)

		deallocate(syTvals)
		deallocate(noSys)

		write(*,*) "Entered solution step"

		call bicgstab(sySt,syRowPtr,syCols,sySrc,100000,revals,iter)
		write(*,'(a,i5,2x,a)') "This program took: ",iter,			&
		"iterations to converge."

		open(resfilenum,file=resfile)

		do i=1,nNodes
			write(resfilenum,'(3(f9.4,2x),f9.4)')noVerts(i,1:3), 	&
			revals(i)
		end do

		close(resfilenum)

		deallocate(noVerts)

	end subroutine fem

    subroutine getmeshdata(meshdetails,vertices,connectivity,		&
	domainelements,surfacenames,surfacefaces)
        integer,parameter :: unitnumber = 111
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

	subroutine setupfinalequations(noSys,sySrc,syTvals)
		integer :: i
		real(8) :: sySrc(:),syTvals(:)
		type(noderow):: noSys(:)

		do i=1,size(noSys,1)
			if(syTvals(i) .ne. 0.0d0) then
				if(allocated(noSys(i)%col)) then
					deallocate(noSys(i)%val)
					deallocate(noSys(i)%col)
				end if
				allocate(noSys(i)%col(1))
				allocate(noSys(i)%val(1))
				noSys(i)%col(1) = i
				noSys(i)%val(1) = 1.0d0
				sySrc(i) = syTvals(i)
			end if
		end do
	end subroutine setupfinalequations

end module htfem
