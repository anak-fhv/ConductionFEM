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
!		integer,parameter :: resfilenum=101, outfilenum=102,		&	! Files to store results and other output
!							 emfilenum=105,nbins=100, bindim=3,		&
!							 srcfilenum=106
		integer :: nNodes,nElems,nDoms,nSurfs,elDom,fcBytype,i,j,k,	&	! Prefixes: n=>number, el=>element, fc=>face
				   typMat,emDom,emFc,iter,meshVals(7),elNodes(4),	&	! by=>boundary, sf=>surface, gn=>generation, no=>node
				   elByfaces(4)
		integer,allocatable :: doElems(:),byCs(:),stRowPtr(:),		&	! do=>domain, sy=>global system
							   stCols(:),cpRowPtr(:),cpCols(:),		&
							   connTab(:,:),sfElems(:,:)
		real(8),parameter :: kDefault = 1.d0
		real(8) :: elVol,tAmbient,gnVal,elK,tc,dlow,dhigh,qBHigh,	&
				   qBLow,abCo,elVolEm,elSurfEm,smFact,byTemp(4),	&
				   bySrc(4),gnSrc(4),elVerts(4,3),elSpfns(4,4),		&
				   elSt(4,4),bySt(4,4),elCp(4,4),noSrcRad
		real(8),allocatable :: domKs(:),sfVals(:),sySt(:),sySrc(:), &
							   syTvals(:),noVerts(:,:),reVals(:),	&
							   vF(:),domRhos(:),domCs(:),syCp(:),	&
							   initGuess(:),rescheck(:),smres(:)
		character(*),parameter :: datdir = "../data/",				&
								  pDatFile = "problemdata.dat",		&
								  resdir = "../results/",			&
								  objdir = "../obj/",				&
								  emfile = objdir//"emissions.out",	&
								  srcfile = objdir//"radsource.out"
		character(72) :: meshFile,byFile,resFile,resVtk,oldResFile,	&
		transResPre
		character(16),allocatable :: sfFcname(:)
		logical,parameter :: gnDefault = .false.,trDefault = .false.
		logical :: trUser,coupled,gnUser,useRK,doSmoothing
		type(noderow) :: noElemPart(4),cpElemPart(4)
		type(noderow),allocatable :: stNo(:),cpNo(:)
		type(elementbins),allocatable :: elbins(:)
		type(surfaceData),allocatable :: surfaces(:)

		call getproblemdata(datdir,pDatFile,meshFile,byFile,		&			! Make an arrangement for material data to be read here: Make the boundary reader exclusive to boundaries.
		resFile,resVtk,trUser,coupled)

		write(*,'(a)') "Now reading mesh file..."
		meshFile = datdir//trim(adjustl(meshFile))

		call getmeshdata(meshFile,meshVals,noVerts,connTab,doElems,	&
		sfFcname,sfElems,surfaces)

		write(*,'(a)') "Meshdetails received."

		nNodes = meshVals(1)
		nElems = meshVals(2)
		nDoms  = meshVals(6)
		nSurfs = meshVals(7)

!		dlow = minval(noverts(:,bindim),1)
!		dhigh = maxval(noverts(:,bindim),1)

		allocate(sySrc(nNodes))
		allocate(syTvals(nNodes))
		allocate(stNo(nNodes))
!		allocate(elbins(nbins))

		if(nDoms .gt. 1) then
			allocate(vF(nDoms))
			vF = 0.d0
		end if

		sySrc = 0.d0
		syTvals = 0.d0

		byFile = datdir//trim(adjustl(byFile))
		call readboundaryconditions(byFile,meshVals,byCs,domKs,		&
		domRhos,domCs,sfVals,tAmbient,gnUser,gnVal)

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

		end do

		write(*,'(a)') "Assembly completed."

		if(trUser) then

! 			Empty call created for later implementation
!			call readinitialvalues(syTvals)

			syTvals = 100.d0
			call collapse_noderows(stNo,sySt,stCols,stRowPtr)
			deallocate(stNo)
			call collapse_noderows(cpNo,syCp,cpCols,cpRowPtr)
			deallocate(cpNo)
			useRK = .true.
			transResPre = resdir//"res_iter_"
			call transientsolve(sySt,stRowPtr,stCols,syCp,cpRowPtr,	&
			cpCols,sySrc,useRK,syTvals,noVerts,connTab,nDoms,		&
			doElems,transResPre)

			return

		end if


		oldResFile = resdir//"results0.out"
		abCo = 1.d0
		call ansourcesurf(trim(adjustl(oldResFile)),abCo,noVerts,	&
		connTab,surfaces(11),sySrc)

		call writepresetupdata(stNo,sySrc)

		call setupfinalequations(stNo,sySrc,syTvals)

!		coupled = .false.
		call collapse_noderows(stNo,sySt,stCols,stRowPtr)
		deallocate(stNo)

		write(*,*)""
		write(*,'(a)') "Entered solution step..."

		call getInitialGuess(syTvals,noVerts,initGuess)
		initGuess = 0.d0

		call bicgstab(sySt,stRowPtr,stCols,sySrc,100000,initGuess,	&
		reVals,iter)

		write(*,'(a)') "Solution completed."
		write(*,'(a,i5,2x,a)')"This program took: ",iter,"iterations&
		& to converge."

		doSmoothing = .true.
		if(doSmoothing) then
!			oldResFile = trim(adjustl(datdir))//oldResFile
			smFact = 0.99d0
			call expsmooth(reVals,smFact,trim(adjustl(oldResFile)),	&
			smRes)
			reVals = smRes
		end if

		call writeresults(noVerts,connTab,nDoms,doElems,reVals,		&
		resdir,resFile,resVtk)

		call writesurfaceemission(noVerts,connTab,reVals,surfaces(11))

		if(coupled) then
			call getelememissions(typMat,emDom,doElems,abCo,reVals,	&
			connTab,sfElems)
		end if

	end subroutine fem

!--------------------------------------------------------------------
!	Subroutine to get essential problem data
!--------------------------------------------------------------------

	subroutine getproblemdata(datdir,pDatFile,meshFile,byFile,		&
	resFile,resVtk,trUser,coupled)
		integer,parameter :: fno = 101
		character(*) :: datdir,pDatFile
		character(72) :: meshFile,byFile,resFile,resVtk
		logical :: trUser,coupled

		open(fno,file=datdir//pDatFile)
		read(fno,*)
		read(fno,*) meshFile
		read(fno,*)
		read(fno,*)
		read(fno,*) byFile
		read(fno,*)
		read(fno,*)
		read(fno,*) resFile
		read(fno,*)
		read(fno,*)
		read(fno,*) resVtk
		read(fno,*)
		read(fno,*)
		read(fno,*) trUser
		read(fno,*)
		read(fno,*)
		read(fno,*) coupled
		close(fno)
	end subroutine getproblemdata

!--------------------------------------------------------------------
!	End subroutine
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Smooth results (exponential smoothing)
!--------------------------------------------------------------------

	subroutine expsmooth(reVals,smFact,oldResFile,smRes)
		integer,parameter :: oldResFno = 105
		integer :: i,nNodes
		real(8) :: smFact,temp(4),reVals(:)
		real(8),allocatable :: smRes(:),oldReVals(:)
		character(*) :: oldResFile

		write(*,*) "Oldresfile: ",oldResFile
		nNodes = size(reVals,1)
		allocate(smRes(nNodes))
		allocate(oldReVals(nNodes))

		open(oldResFno,file=oldResFile)
		do i=1,nNodes
			read(oldResFno,*) temp
			oldReVals(i) = temp(4)
		end do
		close(oldResFno)

		smRes = smFact*oldReVals + (1.d0-smFact)*reVals

		write(*,*) "Compare: "
		write(*,*) oldReVals(5), reVals(5), smRes(5)
	end subroutine expsmooth

!--------------------------------------------------------------------
!	End subroutine
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Verify against PARDISO
!--------------------------------------------------------------------

	subroutine rescheck(noVerts,sySt,sySrc,stRowPtr,stCols,resdir,	&
	resFile,resmkl)
		integer :: i,nNodes,stRowPtr(:),stCols(:)
		integer,parameter :: resChFileNum = 104
		real(8) :: sySt(:),sySrc(:),noVerts(:,:)
		real(8),allocatable :: resmkl(:)
		character(*),parameter :: rFP = "pardiso_"
		character(*) :: resdir,resFile

		nNodes = size(sySrc,1)
		allocate(resmkl(nNodes))
		call solvesystem(sySt,sySrc,stRowPtr,stCols,resmkl)

		resFile = trim(adjustl(resdir))//rFP//trim(adjustl(resFile))
		open(resChFileNum,file=resFile)

		do i=1,nNodes
			write(197,'(4e16.8,2x)')noVerts(i,1:3),resmkl(i)
		end do

		close(197)		
	end subroutine rescheck

!--------------------------------------------------------------------
!	End verification
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Result writer
!--------------------------------------------------------------------

	subroutine writeresults(noVerts,connTab,nDoms,doElems,reVals,	&
	resdir,resFile,resVtk)
		integer,parameter :: resFileNum = 103
		integer :: i,nDoms,nNodes,doElems(:),connTab(:,:)
		real(8) :: reVals(:),noVerts(:,:)
		character(*) :: resdir,resFile,resVtk

		resFile = resdir//trim(adjustl(resFile))
		resVtk = resdir//trim(adjustl(resVtk))
		nNodes = size(reVals,1)

		open(resFileNum,file=resfile)
		do i=1,nNodes
			write(resfilenum,'(4e16.8,2x)')noVerts(i,1:3),revals(i)
		end do
		close(resFileNum)

		call writeresultsvtk(noVerts,connTab,nDoms,doElems,reVals,	&
		resVtk)
		
	end subroutine writeresults

!--------------------------------------------------------------------
!	End result writer
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Logger for stiffness and source terms
!--------------------------------------------------------------------

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

!--------------------------------------------------------------------
!	End logger
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Analytical source term for 2-domain, middle plate case
!--------------------------------------------------------------------

	subroutine ansourcesurf(oldresfile,aC,noVerts,connTab,emSurf,	&
	sySrc)
		integer,parameter :: oldresfno = 106
		integer :: i,j,k,nEmFc,emEl,emElFc,nNodes,ltc,utc,bFcNo(3),	&
		connTab(:,:)
		real(8) :: aC,T1,T2,T3,Tcent,fcA,emVal,Tl,Th,recVal,		&
		sfSrc(3),cent(3),temp(4),noVerts(:,:)
		real,allocatable :: reVals(:)
		real(8),parameter :: sigb = 5.670373e-8, kel = 273.15d0,	&
		Tb1 = 373.15d0, Tb2 = 283.15d0
		real(8),intent(inout) :: sySrc(:)
		character(*) :: oldresfile
		type(surfaceData) :: emSurf


		nNodes = size(sySrc,1)
		allocate(reVals(nNodes))
		open(oldresfno,file=oldresfile)
		do i=1,nNodes
			read(oldresfno,*) temp
			reVals(i) = temp(4)
		end do
		close(oldresfno)

		nEmFc = size(emSurf%elNum,1)
		aC = 1.d0

		ltc = 0
		utc = 0
		do i=1,nEmFc
			emEl = emSurf%elNum(i)
			emElFc = emSurf%fcNum(i)
			call bfacenodes(emElFc,bFcNo)
			bFcNo = connTab(emEl,bFcNo)
			fcA = facearea(noVerts(bFcNo,:))
			T1 = reVals(bFcNo(1))
			T2 = reVals(bFcNo(2))
			T3 = reVals(bFcNo(3))
			Tcent = (T1+T2+T3)/3.d0
			if(Tcent .lt. 0.d0) then
				Tcent = 0.d0
			end if
			emVal = fcA*sigB*aC*(Tcent**4.d0)
			cent = sum(noVerts(bFcNo,3))/3.d0
			if(cent(3).lt. 0.d0) then
				ltc = ltc+1
				Tl = Tl + Tcent
				recVal = aC*sigb*(Tb1**4.d0)*fcA
				sfSrc = (recVal-emVal)/3.d0
				sySrc(bFcNo) = sySrc(bFcNo) + sfSrc
			else
				utc = utc+1
				Th = Th + (T1+T2+T3)/3.d0
				recVal = aC*sigb*(Tb2**4.d0)*fcA
				sfSrc = (recVal-emVal)/3.d0
				sySrc(bFcNo) = sySrc(bFcNo) + sfSrc
			end if
		end do

	end subroutine ansourcesurf

!--------------------------------------------------------------------
!	End analytical source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Element emissions based on case
!--------------------------------------------------------------------

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
								call getsurfaceemission(abCo,		&
								reVals(elNo),emFc,elSurfEm)
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

!--------------------------------------------------------------------
!	End element emissions
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!	Volumetric emission and surface source - analytical source calc
!--------------------------------------------------------------------

	subroutine ansourcevol(oldresfile,aC,noVerts,connTab,bSurfs,	&
	sySrc)
		integer,parameter :: oldresfno = 107
		integer :: i,j,k,nEmFc,nElems,emEl,emElFc,nNodes,ltc,utc,	&
		bFcNo(3),connTab(:,:)
		real(8) :: aC,T1,T2,T3,T4,Tcent,fcA,emVal,emPerVol,Tl,Th,	&
		recVal,sfSrc(3),cent(3),temp(4),noVerts(:,:)
		real,allocatable :: reVals(:)
		real(8),parameter :: sigb = 5.670373e-8, kel = 273.15d0,	&
		Tb1 = 373.15d0, Tb2 = 283.15d0
		real(8),intent(inout) :: sySrc(:)
		character(*) :: oldresfile
		type(surfaceData) :: bSurfs(:)

		nEmSurf = size(bSurfs,1)
		nElems = size(connTab,1)
		emVal = 0.d0
		do i=1,nEmSurf
			nEmFc = size(bSurfs(i)%elNum,1)
			do j=1,nEmFc
				emEl = emSurf%elNum(i)
				emElFc = emSurf%fcNum(i)
				call bfacenodes(emElFc,bFcNo)
				bFcNo = connTab(emEl,bFcNo)
				fcA = facearea(noVerts(bFcNo,:))
				T1 = reVals(bFcNo(1))
				T2 = reVals(bFcNo(2))
				T3 = reVals(bFcNo(3))
				Tcent = (T1+T2+T3)/3.d0
				emVal = emVal + fcA*sigB*aC*(Tcent**4.d0)
			end do
		end do

		do i=1,nElems
			elNodes = connTab(i,:)
			elVerts = noVerts(elNodes,:)
						
		end do

	subroutine ansourcevol

!--------------------------------------------------------------------
!	End volumetric sources
!--------------------------------------------------------------------

    subroutine getmeshdata(meshFile,meshdetails,vertices,			&
	connectivity,domainelements,surfacenames,surfacefaces,surfaces)
        integer,parameter :: unitnumber = 102
		integer :: meshdetails(7)
        integer,allocatable :: domainelements(:),connectivity(:,:),	&
		surfacefaces(:,:)
		character(*) :: meshFile
        character(len=16),allocatable :: surfacenames(:)
        real(8),allocatable :: boundaryvalues(:),vertices(:,:)
		type(surfaceData),allocatable :: surfaces(:)

        call openmeshfile(unitnumber, meshFile)
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
