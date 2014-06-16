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
		integer,parameter :: resfilenum=101, outfilenum=102,		&	! Files to store results and other output
							 emfilenum=105,nbins=100, bindim=3
		integer :: nNodes,nElems,nDoms,nSurfs,elDom,fcBytype,i,j,k,	&	! Prefixes: n=>number, el=>element, fc=>face
				   typMat,emDom,emFc,iter,meshVals(7),elNodes(4),	&	! by=>boundary, sf=>surface, gn=>generation, no=>node
				   elByfaces(4)
		integer,allocatable :: doElems(:),byCs(:),stRowPtr(:),		&	! do=>domain, sy=>global system
							   stCols(:),cpRowPtr(:),cpCols(:),		&
							   connTab(:,:),sfElems(:,:)
		real(8),parameter :: kDefault = 1.d0
		real(8) :: elVol,tAmbient,gnVal,elK,tc,dlow,dhigh,qBHigh,	&
				   qBLow,absCoeff,elVolEm,elSurfEm,byTemp(4),		&
				   bySrc(4),gnSrc(4),elVerts(4,3),elSpfns(4,4),		&
				   elSt(4,4),bySt(4,4),elCp(4,4)
		real(8),allocatable :: domKs(:),sfVals(:),sySt(:),sySrc(:), &
							   syTvals(:),noVerts(:,:),reVals(:),	&
							   vF(:),domRhos(:),domCs(:),syCp(:),	&
							   initGuess(:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"results.out",	&
								  outfile = objdir//"outputs.out",	&
								  emfile = objdir//"emissions.out"
		character(16),allocatable :: sfFcname(:)
		logical,parameter :: gnDefault = .false.,trDefault = .false.
		logical ::	gnUser,trUser,useRK,radUser
		type(noderow) :: noElemPart(4),cpElemPart(4)
		type(noderow),allocatable :: stNo(:),cpNo(:)
		type(elementbins),allocatable :: elbins(:)
		type(surfaceData),allocatable :: surfaces(:)

		write(*,'(a)') "Now reading mesh file..."

		call getmeshdata(meshVals,noVerts,connTab,doElems,			&
		sfFcname,sfElems,surfaces)

		write(*,'(a)') "Meshdetails received."

		nNodes = meshVals(1)
		nElems = meshVals(2)
		nDoms  = meshVals(6)
		nSurfs = meshVals(7)

		dlow = minval(noverts(:,bindim),1)
		dhigh = maxval(noverts(:,bindim),1)

		allocate(sySrc(nNodes))
		allocate(syTvals(nNodes))
		allocate(stNo(nNodes))
		allocate(elbins(nbins))

		if(nDoms .gt. 1) then
			allocate(vF(nDoms))
			vF = 0.d0
		end if

		sySrc = 0.d0
		syTvals = 0.d0

		call readboundaryconditions(meshVals,byCs,domKs,domRhos,	&
		domCs,sfVals,tAmbient,trUser,gnUser)

!		trUser = .false.
		if(trUser) then
			allocate(cpNo(nNodes))
		end if

		call summarisesystem(meshVals,byCs,sfVals,sfFcname,trUser,	&
		gnUser)

		write(*,'(a)') "Beginning assembly..."

		do i=1,nElems

!	Get element data from mesh values
			elNodes = connTab(i,:)
			elDom = doElems(i)
			elVerts = noVerts(elNodes,:)
			elK = domKs(elDom)

!	Call stiffness functions
			call shapefunctions(elVerts,elVol,elSpfns)
			call elementstiffness(elSpfns,elVol,elK,elSt)

!			call binelement(elbins,bindim,dlow,dhigh,i,elVerts)

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

			if(trUser) then
				cpElemPart = cpNo(elNodes)
				call assemblecapacitance(cpElemPart,elNodes,elCp)
				cpNo(elNodes) = cpElemPart
			end if
		end do

		write(*,'(a)') "Assembly completed."

		if(trUser) then
! Empty call created for later implementation
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

		if(trUser) then
			deallocate(cpNo)
		end if

		write(*,*)""
		write(*,'(a)') "Entered solution step..."

		if(trUser) then
!			useRK = .true.
			call transientsolve(sySt,stRowPtr,stCols,syCp,cpRowPtr,	&
			cpCols,sySrc,useRK,syTvals,noVerts,connTab,nDoms,		&
			doElems)
		else
			
			call getInitialGuess(syTvals,noVerts,initGuess)
			call bicgstab(sySt,stRowPtr,stCols,sySrc,100000,		&
			initGuess,reVals,iter)

			write(*,'(a)') "Solution completed."
			write(*,'(a,i5,2x,a)') "This program took: ",iter,		&
			"iterations to converge."

!--------------------------------------------------------------------
!	Function to check the effect of face area on the emission value
!--------------------------------------------------------------------
!			call checkemissiondifference(noVerts,connTab,sfElems,	&
!			13,reVals)


!--------------------------------------------------------------------
!	Function to write facewise emission values for a given surface
!--------------------------------------------------------------------

			call writesurfaceemission(noVerts,connTab,reVals,		&
			surfaces(13))

			open(resfilenum,file=resfile)

			do i=1,nNodes
				write(resfilenum,'(3(f9.4,2x),f9.4)')noVerts(i,1:3),&
				revals(i)
			end do

!			write(resfilenum,*)

!			call getflowrates(noVerts,connTab,doElems,domKs,sfElems,&
!			reVals,(/1/),(/5/),3,qBLow,qBHigh)

!			if(nDoms .eq. 2) then
!				write(resfilenum,*)"Sample porosity:",vF(1)/(sum(vF))
!			end if

!			write(resfilenum,*) "Fluxes:"
!			write(resfilenum,*) "Boundary low:", qBLow
!			write(resfilenum,*) "Boundary high:", qBHigh

			close(resfilenum)

!--------------------------------------------------------------------
!	Function to get element emissions
!--------------------------------------------------------------------

!			if(radUser) then
!				open(emfilenum,file=emfile)
!				if(typMat == 1) then
!					do i=1,nElems
!						if(doElems(i)==emDom) then
!							call getelementvolumeemission(absCoeff,	&
!							reVals(connTab(i,:)),elVolEm)
!						end if
!						write(emfilenum,*) i,elVolEm
!					end do
!				elseif(typMat == 2) then
!					do i=1,nElems
!						if(doElems(i) == emDom) then
!							if(any(sfElems(i,:) == emSurf)) then
!								emFc = sfElems(minloc(sfElems,		&
!								sfElems==emSurf))
!								call getsurfaceemission(absCoeff,	&
!								reVals(connTab(i,:)),emFc,elSurfEm)
!							end if
!						end if
!						write(emfilenum,*) i,elSurfEm
!					end do
!				else
!					write(*,'(a,1x,i3,1x,a)') "The material type &
!					&specified is: ", typMat, "an unknown. Please &
!					&update the data file to show a correct &
!					&material type."
!				end if
!				close(emfilenum)
!			end if

			call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
			reVals)
		end if

	end subroutine fem

    subroutine getmeshdata(meshdetails,vertices,connectivity,		&
	domainelements,surfacenames,surfacefaces,surfaces)
        integer,parameter :: unitnumber = 103
		integer :: meshdetails(7)
        integer,allocatable :: domainelements(:),connectivity(:,:),	&
		surfacefaces(:,:)
        character(len=16),allocatable :: surfacenames(:)
        real(8),allocatable :: boundaryvalues(:),vertices(:,:)
		type(surfaceData),allocatable :: surfaces(:)

        call openmeshfile(unitnumber, 'a.msh')
        call readmeshdetails(unitnumber,meshdetails)
        call readmeshvertices(unitnumber,meshdetails, vertices)
        call readmeshconnectivity(unitnumber,meshdetails, 			&
		connectivity)
        call readmeshdomains(unitnumber,meshdetails,domainelements)
        call readmeshsurfaces(unitnumber,meshdetails,surfacefaces,	&
		surfacenames,surfaces)
        call closemeshfile(unitnumber)
    end subroutine getmeshdata

	subroutine summarisesystem(meshVals,byCs,sfVals,sfFcname,		&
	trUser,gnUser)
		integer :: i,meshVals(:),byCs(:)
		real(8),parameter :: small=1e-8
		real(8) :: sfVals(:)
		logical :: trUser,gnUser
		character(*) :: sfFcname(:)

		write(*,*)""
		write(*,'(a)') "Short system summary: "
		write(*,'(4x,a,1x,i4)') "Domains: ", meshVals(6)
		write(*,'(4x,a,1x,i8)') "Number of nodes: ", meshVals(1)
		write(*,'(4x,a,1x,i8)') "Number of elements: ", meshVals(2)
		write(*,*)""
		write(*,'(4x,a)') "System boundaries, not including &
		&domain interfaces: "
		write(*,*)""

		do i=1,meshVals(7)
			
			if(byCs(i) == 1) then
				write(*,'(4x,a,1x,a)') trim(sfFcname(i)), &
				&": Temperature"
			elseif(byCs(i) == 2) then
				if(sfVals(i) < small) then
					write(*,'(4x,a,1x,a)') trim(sfFcname(i)), &
					&": Adiabatic"
				else
					write(*,'(4x,a,1x,a)') trim(sfFcname(i)), &
					&": Flux"
				end if
			elseif(byCs(i) == 3) then
				write(*,'(4x,a,1x,a)') trim(sfFcname(i)), &
				&": Convective"
			elseif(byCs(i) == 4) then
				continue
			else
				write(*,'(4x,a)') "For ", trim(sfFcname(i))
				write(*,'(4x,a)')"Unrecognised boundary condition."
				write(*,'(4x,a)')"Update datafile.dat inputs."
				write(*,'(4x,a)')"Now quitting program execution..."
				call exit(1)
			end if
		end do

		if(gnUser) then
			write(*,'(4x,a)') "User specified volume generation&
			& exists."
		end if

		if(trUser) then
			write(*,'(4x,a)') "Transient solution requested."
		end if

		write(*,*)""
	end subroutine summarisesystem

	subroutine binelement(elbins,bindim,dlow,dhigh,elno,ec)
		integer :: bindim,elno
		real(8) :: dlow,dhigh,cent(3),ec(4,3)
		type(elementbins) :: elbins(:)

		call elementcentroid(ec,cent)
		call addtoelementbins(elno,cent,bindim,dlow,dhigh,elbins)
	end subroutine binelement

	subroutine addtoglobaltemperature(Tvals,elnodes,btemp)
		real(8),intent(inout) :: Tvals(:)
		real(8) :: btemp(4)
		integer :: elnodes(4)

		if(any(btemp.ne.0.0d0)) then
			Tvals(elnodes) = btemp
		end if
	end subroutine addtoglobaltemperature

	subroutine addtoglobalforce(gF,elnodes,elf)
		real(8),intent(inout) :: gF(:)
		real(8) :: elf(4)
		integer :: i,elnodes(4)

		do i=1,4
			gF(elnodes(i)) = gF(elnodes(i)) + elf(i)
		end do
	end subroutine addtoglobalforce

	subroutine uniformgeneration(gval,elvol,gcontrib)
		real(8) :: gval,elvol,gcontrib(4)

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

	subroutine getInitialGuess(syTvals,noVerts,initGuess)
		integer :: i,j,n,pos
		real(8) :: vMax,vMin,Tmax,Tmin,grad,syTvals(:),noVerts(:,:)
		real(8),allocatable :: initGuess(:),tVerts(:)

		n = size(noVerts,1)
		allocate(tVerts(n))
		allocate(initGuess(n))

		tVerts = noVerts(:,3)
		vMax = maxval(tVerts)
		vMin = minVal(tVerts)

		write(*,*) "maxloc: ", maxloc(tVerts)
		write(*,*) "minloc: ", minloc(tVerts)
		Tmax = syTvals(maxloc(tVerts,1))
		Tmin = syTvals(minloc(tVerts,1))
		write(*,*) "maxval: ", Tmax
		write(*,*) "minval: ", Tmin
		if(Tmax .ne. Tmin) then
			grad = (Tmax-Tmin)/(vMax-vMin)
			do i=1,n
				initGuess(i)=Tmin+(tVerts(i)-vMin)*grad
			end do
		else
			initGuess = 0.d0
		end if

	end subroutine getInitialGuess

end module htfem
