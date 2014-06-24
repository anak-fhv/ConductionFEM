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
							 emfilenum=105,nbins=100, bindim=3,		&
							 srcfilenum=106
		integer :: nNodes,nElems,nDoms,nSurfs,elDom,fcBytype,i,j,k,	&	! Prefixes: n=>number, el=>element, fc=>face
				   typMat,emDom,emFc,iter,meshVals(7),elNodes(4),	&	! by=>boundary, sf=>surface, gn=>generation, no=>node
				   elByfaces(4)
		integer,allocatable :: doElems(:),byCs(:),stRowPtr(:),		&	! do=>domain, sy=>global system
							   stCols(:),cpRowPtr(:),cpCols(:),		&
							   connTab(:,:),sfElems(:,:)
		real(8),parameter :: kDefault = 1.d0
		real(8) :: elVol,tAmbient,gnVal,elK,tc,dlow,dhigh,qBHigh,	&
				   qBLow,abCo,elVolEm,elSurfEm,byTemp(4),bySrc(4),	&
				   gnSrc(4),elVerts(4,3),elSpfns(4,4),elSt(4,4),	&
				   bySt(4,4),elCp(4,4),noSrcRad
		real(8),allocatable :: domKs(:),sfVals(:),sySt(:),sySrc(:), &
							   syTvals(:),noVerts(:,:),reVals(:),	&
							   vF(:),domRhos(:),domCs(:),syCp(:),	&
							   initGuess(:),rescheck(:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"results.out",	&
								  outfile = objdir//"outputs.out",	&
								  emfile = objdir//"emissions.out",	&
								  srcfile = objdir//"radsource.out"
		character(16),allocatable :: sfFcname(:)
		logical,parameter :: gnDefault = .false.,trDefault = .false.
		logical ::	gnUser,trUser,useRK,radUser=.false.
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

!	Random source terms in the elements

!			tval = ?
!			call randomradiationsource(abCo,emVal,tVal,elNodes,		&
!			sySrc)

		end do

		write(*,'(a)') "Assembly completed."

		if(trUser) then
! Empty call created for later implementation
!			call readinitialvalues(syTvals)
			syTvals = 100.d0
		else
!!			radUser = .true.
!			if(radUser) then
!				open(srcfilenum,file=srcfile)
!				do i=1,nNodes
!					read(srcfilenum,*) noSrcRad
!					sySrc(i) = sySrc(i)+noSrcRad
!				end do
!				close(srcfilenum)
!			end if

!			allocate(rescheck(nNodes))
!			open(2468,file=objdir//"resvals.out")
!			do i=1,nNodes
!				read(2468,*) byTemp
!				rescheck(i) = byTemp(4)
!			end do
!			close(2468)
!			call ansource(noVerts,connTab,rescheck,surfaces(11),sySrc)

!			deallocate(rescheck)

!			open(1357,file=objdir//"tempsource.out")
!			write(1357,'(e20.8)') sySrc
!			close(1357)

!			call writepresetupdata(stNo,sySrc)

			call setupfinalequations(stNo,sySrc,syTvals)
		end if

		radUser = .false.
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
			initGuess = 0.d0

			call bicgstab(sySt,stRowPtr,stCols,sySrc,100000,		&
			initGuess,reVals,iter)

			open(resfilenum,file=resfile)

			do i=1,nNodes
				write(resfilenum,'(3(f9.4,2x),e20.8)')noVerts(i,1:3),&
				revals(i)
			end do

			close(resfilenum)

!			allocate(rescheck(nNodes))
!			call solvesystem(sySt,sySrc,stRowPtr,stCols,rescheck)

!			open(197,file=objdir//"rescomp.out")

!			do i=1,nNodes
!				write(197,'(3(f9.4,2x),e20.8)')noVerts(i,1:3),&
!				rescheck(i)
!			end do

!			close(197)

			write(*,'(a)') "Solution completed."
			write(*,'(a,i5,2x,a)') "This program took: ",iter,		&
			"iterations to converge."

!--------------------------------------------------------------------
!	Function to write facewise emission values for a given surface
!--------------------------------------------------------------------

			call writesurfaceemission(noVerts,connTab,reVals,		&
			surfaces(11))

			open(resfilenum,file=resfile)

			do i=1,nNodes
				write(resfilenum,'(3(f9.4,2x),e20.8)')noVerts(i,1:3),&
				revals(i)
			end do

			close(resfilenum)

!--------------------------------------------------------------------
!	Function to get element emissions
!--------------------------------------------------------------------
			if(radUser) then
				call getelememissions(typMat,emDom,doElems,abCo,	&
				reVals,connTab,sfElems)
			end if
!--------------------------------------------------------------------
!	Function to get element emissions
!--------------------------------------------------------------------

			call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
			reVals)
		end if

	end subroutine fem

	subroutine writepresetupdata(stNo,sySrc)
		integer,parameter :: valfilno = 195,srcfilno = 153
		integer :: i,j,k,nNodes,npercol
		real(8) :: sySrc(:)
		character(*),parameter :: objdir = "../obj/",				&
								  valfil = objdir//"stvals.out",	&
								  srcfil = objdir//"scvals.out"
		type(noderow) :: stNo(:)

		nNodes = size(stNo,1)
		open(valfilno,file=valfil)
		open(srcfilno,file=srcfil)
		do i=1,nNodes
			write(srcfilno,*) sySrc(i)
			write(valfilno,*) "row: ",i
			npercol = size(stNo(i)%col,1)
			do j=1,npercol
				write(valfilno,'(i8,2x,e20.8)') stNo(i)%col(j),stNo(i)%val(j)
			end do
		end do
	end subroutine writepresetupdata

	subroutine ansource(noVerts,connTab,reVals,emSurf,sySrc)
		integer :: i,j,k,nEmFc,emEl,emElFc,bFcNo(3),connTab(:,:)
		real(8) :: aC,T1,T2,T3,fcA,emVal,recVal,reVals(:),			&
		noVerts(:,:),sfSrc(3),cent(3)
		real(8),parameter :: sigb = 5.670373e-8, kel = 273.15d0,	&
		Th = 373.15d0, Tl = 283.15d0
		real(8),intent(inout) :: sySrc(:)
		type(surfaceData) :: emSurf

		nEmFc = size(emSurf%elNum,1)
		aC = 1.d0
		do i=1,nEmFc
			emEl = emSurf%elNum(i)
			emElFc = emSurf%fcNum(i)
			call bfacenodes(emElFc,bFcNo)
			bFcNo = connTab(emEl,bFcNo)
			fcA = facearea(noVerts(bFcNo,:))
			T1 = reVals(bFcNo(1))		! + kel
			T2 = reVals(bFcNo(2)) 		! + kel
			T3 = reVals(bFcNo(3)) 		! + kel
			emVal = fcA*sigB*aC*((T1+T2+T3)/3.d0)**4
			cent = sum(noVerts(bFcNo,:))/3.d0
			if(cent(3).lt.0.d0) then
				recVal = aC*sigb*(Th**4.d0)*fcA
				sfSrc = (recVal-emVal)/3.d0
				sySrc(bFcNo) = sySrc(bFcNo) + sfSrc
			else
				recVal = aC*sigb*(Tl**4.d0)*fcA
				sfSrc = (recVal-emVal)/3.d0
				sySrc(bFcNo) = sySrc(bFcNo) + sfSrc
			end if
		end do

	end subroutine ansource

	subroutine getelememissions(typMat,emDom,doElems,abCo,reVals,	&
	connTab,sfElems)
		integer,parameter :: emFNo = 147
		integer :: i,j,k,typMat,emDom,nElems,emFc,emSurf,elNo(4),	&
		doElems(:),connTab(:,:),sfElems(:,:)
		real(8) :: abCo,elVolEm,elSurfEm,reVals(:)
		character(*),parameter :: emFile="../obj/faceEmissions.out"

		
		nElems = size(doElems,1)
		open(emFNo,file=emFile)
		if(typMat == 1) then
			do i=1,nElems
				elNo = connTab(i,:)
				if(doElems(i)==emDom) then
					call getelementvolumeemission(abCo,reVals(elNo),&
					elVolEm)
				end if
				write(emFNo,*) i,elVolEm
			end do
		elseif(typMat == 2) then
			do i=1,nElems
				elNo = connTab(i,:)
				if(doElems(i) == emDom) then
					if(any(sfElems(i,:) == emSurf)) then
						do j=1,4
							if(sfElems(i,j) == emSurf) then
								emFc = j
								call getsurfaceemission(abCo,reVals(elNo),emFc,elSurfEm)
								write(emFNo,*) i,elSurfEm
							end if
						end do
					end if
				end if

			end do
		else
			write(*,'(a,1x,i3,1x,a)') "The material type &
			&specified is: ", typMat, "an unknown. Please &
			&update the data file to show a correct &
			&material type."
		end if

		close(emFNo)
	end subroutine getelememissions

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

	subroutine randomradiationsource(threshold,emVal,tVal,elNodes,	&
	sySrc)
		integer :: elNodes(4)
		real(8),parameter :: sigb=5.670373e-8
		real(8) :: threshold,check,emVal,tVal,radSource(4)
		real(8),intent(inout) :: sySrc(:)

		call random_number(check)
		if(check.ge.threshold) then
			radSource = sigb*emVal*(tVal**4.d0)
			sySrc(elNodes) = sySrc(elNodes)+radSource
		end if
	end subroutine randomradiationsource

end module htfem
