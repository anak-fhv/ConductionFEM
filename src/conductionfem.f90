include 'mkl_pardiso.f90'	! Added for pardiso

module conductionfem

	use readmesh
	use elementcalculations
! Added for omp
	use omp_lib
! End added for omp
	implicit none
	contains
	
	subroutine finitesolution()
		integer,parameter :: resfilnum=999
		integer :: numno,numel,numdo,numsu,		&				! Meshdetails
				i,j,k,eldom,bface,fbtype						! Counters
		integer,dimension(7) :: mdets							! Always the same size
		integer,dimension(4) :: elnodes,bofcs					! Tetrahedral mesh
		integer,dimension(3) :: fcnodes							! Triangular faces
		integer,dimension(:),allocatable :: domel,bcs
		integer,dimension(:,:),allocatable :: conn,sfacs
		double precision :: elvol,bv,Tamb,gval,hval,elk			! From the data file
		double precision,dimension(4) :: btemp,bforce,hta,gcontrib
		double precision,dimension(4,3) :: ec
		double precision,dimension(4,4) :: sp,btdb,bstiff
		double precision,dimension(:),allocatable :: kvals,bvals,gF,Tvals,newF
		double precision,dimension(:,:),allocatable :: ver,gK,newK
		character(len=*),parameter :: outputdir="../obj/",resfil=outputdir//"results.out"
		character(len=16), dimension(:), allocatable :: sfacnam
        logical :: gen
!	Added for CSR
		double precision,dimension(:),allocatable :: acsr,solns
		integer,dimension(:),allocatable :: ia,ja
!	End added for CSR

		call getmeshdata(mdets,ver,conn,domel,sfacnam,sfacs)
		print *, "meshdetails received"
		numno = mdets(1)
		numel = mdets(2)
		numdo = mdets(6)
		numsu = mdets(7)
		allocate(gF(numno))
		allocate(Tvals(numno))
		allocate(gK(numno,numno))
		gF = 0.0d0
		gK = 0.0d0
		Tvals = 0.0d0
        call readboundaryconditions(mdets,bcs,kvals,bvals,Tamb,gen)
		open(resfilnum,file=resfil)

		!$omp parallel &
		!$omp shared (numel,domel,ver,conn,bcs,bvals,kvals,Tamb,gen,sfacs,gval,gK,gF) &
		!$omp private (i,elnodes,eldom,ec,elk,elvol,sp,btdb,bofcs,bstiff,bforce,btemp,gcontrib)
		!$omp do
		do i=1,numel

! Essential part: Stiffness matrix and adding to the global matrix
			elnodes = conn(i,:)
			eldom = domel(i)
			ec = ver(elnodes,:)
			elk = kvals(eldom)
			call shapefunctions(ec,elvol,sp)
			call elementstiffness(sp,elvol,elk,btdb)
			call addtoglobalstiffness(gK,elnodes,btdb)

! Boundary values - computed as a sum per element
			bofcs = sfacs(i,:)
			if(all(bofcs==0)) then
				continue
			else
				call boundaryconditions(ec,bofcs,bcs,bvals,Tamb,bstiff,bforce,btemp)
				call addtoglobalstiffness(gK,elnodes,bstiff)
				call addtoglobaltemperature(Tvals,elnodes,btemp)
				call addtoglobalforce(gF,elnodes,bforce)
			endif

! Accounting for uniform generation, if present within the volume
			if(gen) then
				call uniformgeneration(gval,elvol,gcontrib)
				call addtoglobalforce(gF,elnodes,gcontrib)
			end if
		end do
		!$omp end do
		!$omp end parallel

! Set up the final matrices for the temperature solution at the nodes
		call setupfinalequations(gK,Tvals,gF,newK,newF)

! Trial call to converting the stiffness matrix into a CSR representation
		call makestiffnesscsr(newK,newF,acsr,ia,ja)

! Added code for the use of MKL PARDISO
		call solvesystem(acsr,newF,ia,ja,solns)

! Solve the final equations using the BLAS/LAPACK/MKL libraries
! Temporarily blocked
!		call solvefinalequations(newK,newF)
! End temporarily blocked

! Bookkeeping (writing values of important results)
		open(resfilnum,file=resfil)
		write(resfilnum,'(a)') "Final results: "
		do i=1,numno
			write(resfilnum,'(3(f9.4,2x),f9.4)')ver(i,1:3), solns(i)	! newF(i) changed to solns(i)
		end do
		write(resfilnum,'(a)') "End of final results."
		close(resfilnum)

	end subroutine finitesolution

! Here begin the helper routines

!1. Routine to get mesh data from the readmesh module
    subroutine getmeshdata(meshdetails,vertices,connectivity,domainelements,surfacenames,surfacefaces)
        integer,parameter :: unitnumber = 111
		integer,dimension(7) :: meshdetails
        integer,dimension(:,:),allocatable :: connectivity, surfacefaces
        integer,dimension(:),allocatable :: domainelements
        character(len=16),dimension(:),allocatable :: surfacenames
        double precision,dimension(:,:),allocatable :: vertices
        double precision,dimension(:),allocatable :: boundaryvalues

        call openmeshfile(unitnumber, 'a.msh')
        call readmeshdetails(unitnumber, meshdetails)
        call readmeshvertices(unitnumber, meshdetails, vertices)
        call readmeshconnectivity(unitnumber, meshdetails, connectivity)
        call readmeshdomains(unitnumber, meshdetails, domainelements)
        call readmeshsurfaces(unitnumber,meshdetails,surfacefaces,surfacenames)
        call closemeshfile(unitnumber)
    end subroutine getmeshdata

!2. Routine to read boundary conditions from the data file
    subroutine readboundaryconditions(meshdetails, boundaryconditions, &
    kvals, boundaryvalues, Tambient, generation, generationrate)
        integer,parameter :: datafilenum=222
        character(len=*),parameter :: filename='datafile.dat'
        integer,dimension(7) :: meshdetails
        integer :: numdomains,numboundaries,i,j
        logical :: generation
        integer,dimension(:),allocatable :: boundaryconditions
        double precision :: Tambient
		double precision,optional :: generationrate
        double precision,dimension(:),allocatable :: kvals,boundaryvalues

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

!3. Routine for assembling each element's various stiffness components into the global stiffness
	subroutine addtoglobalstiffness(gK,elnodes,elst)
		double precision,dimension(:,:),intent(inout) :: gK
		double precision,dimension(4,4) :: elst
		integer,dimension(4) :: elnodes
		integer :: i,j,n1,n2

		do i=1,4
			n1 = elnodes(i)
			do j=1,4
				n2 = elnodes(j)
				gK(n1,n2) = gK(n1,n2) + elst(i,j)
			end do
		end do
	end subroutine addtoglobalstiffness

!4. Routine for assembling temperatures
	subroutine addtoglobaltemperature(Tvals,elnodes,btemp)
		double precision,dimension(:),intent(inout) :: Tvals
		double precision,dimension(4) :: btemp
		integer,dimension(4) :: elnodes

		if(any(btemp.ne.0.0d0)) then
			Tvals(elnodes) = btemp
		end if
	end subroutine addtoglobaltemperature

!5. Routine for assembly of the global force vector
	subroutine addtoglobalforce(gF,elnodes,elf)
		double precision,dimension(:),intent(inout) :: gF
		double precision,dimension(4) :: elf
		integer,dimension(4) :: elnodes
		integer :: i

		do i=1,4
			gF(elnodes(i)) = gF(elnodes(i)) + elf(i)
		end do
	end subroutine addtoglobalforce

!6. Giant subroutine to refactor the handling of boundaries from the element counter
	subroutine boundaryconditions(ec,bofcs,bcs,bvals,Tamb,bstiff,bforce,btemp)
		double precision,dimension(4,3) :: ec
		double precision,dimension(4,4) :: bstiff,elht
		double precision,dimension(4) :: bforce,btemp,elq,hta
		double precision,dimension(:),intent(in) :: bvals
		double precision :: Tamb,bv
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

!7. Routine to handle a uniform generation term
	subroutine uniformgeneration(gval,elvol,gcontrib)
		double precision :: gval,elvol
		double precision,dimension(4) :: gcontrib

		gcontrib = gval*(elvol/24.0d0)*(/1,1,1,1/)
	end subroutine uniformgeneration

!8. Routine to manipulate the final global stiffness and global force matrices for solution
	subroutine setupfinalequations(gK,Tvals,gF,newK,newF)
		double precision,dimension(:,:),intent(in) :: gK
		double precision,dimension(:),intent(in) :: Tvals,gF
		double precision,dimension(:,:),allocatable :: newK
		double precision,dimension(:),allocatable :: newF
		double precision,parameter :: small=1e-6
		integer :: i,j,k,knT,unT,n
		
		n=size(gK,1)
		allocate(newK(n,n))
		allocate(newF(n))
		newF = gF
		newK = gK
		knT = 0
		unT = 0

		!$omp parallel &
		!$omp shared (Tvals,knT,unT,newK,newF) &
		!$omp private (i)

		!$omp do
		do i=1,n
			if(abs(Tvals(i)) < small) then
				unT = unT+1
			else
				knT = knT+1
				newF(i) = Tvals(i)
				newK(i,:) = 0.0d0
				newK(i,i) = 1.0d0
			end if
		end do
		!$omp end do
		!$omp end parallel
		print *, "Known and unknown: ", knT, unT
	end subroutine setupfinalequations

!9. Routine to call MKL functions to solve the final system
	subroutine solvefinalequations(K,F)
		double precision,dimension(:,:) :: K
		double precision,dimension(:) :: F
		integer :: m,n,nrhs,lda,ldb,info
		integer,dimension(size(K,1)) :: ipiv
		character(len=1) :: trans="N"

		m = size(K,1)
		n = size(K,1)
		lda = size(K,1)
		ldb = size(K,1)
		nrhs = 1
		call dgetrf(m,n,K,lda,ipiv,info)
		call dgetrs(trans,n,nrhs,K,lda,ipiv,F,ldb,info)
	end subroutine solvefinalequations

!10. Convert the global stiffness into a sparse column representation
	subroutine makestiffnesscsr(K,F,acsr,ia,ja)
		integer,dimension(8) :: job
		integer :: i,j,ct,m,n,lda,info
		integer,dimension(:),allocatable :: ia,ja
		double precision,dimension(:,:) :: K
		double precision,dimension(:) :: F
		double precision,dimension(:),allocatable :: acsr

		m = size(K,1)
		n = size(K,2)
		lda = m

		ct = 0

		!$omp parallel &
		!$omp shared (K,ct,m,n) &
		!$omp private (i,j)
		!omp do
		do i=1,m
			!$omp do
			do j=1,n
				if(K(i,j)/=0.0d0) then
					ct = ct+1
				end if
			end do
			!$omp end do
		end do
		!omp end do
		!$omp end parallel
		print *, "Count: ", ct
		allocate(ja(ct))
		allocate(acsr(ct))
		allocate(ia(m+1))
		job(1:6) = (/0,1,1,2,ct,3/)
		call mkl_ddnscsr(job,m,n,K,lda,acsr,ja,ia,info)
		print *, "info: ", info
	end subroutine makestiffnesscsr

	subroutine solvesystem(Kcsr,F,ia,ja,x)
		use mkl_pardiso
		integer,parameter :: dp = kind(1.0d0)		! internal solver memory pointer 
		type(mkl_pardiso_handle),dimension(:),allocatable  :: pt		! all other variables
		integer :: i,maxfct,mnum,mtype,phase,n,nrhs,error,msglvl,nnz,error1
		integer,dimension(:),allocatable :: iparm
		integer,dimension(:),intent(in) :: ia,ja
		double precision,dimension(:) :: Kcsr,F
		double precision,dimension(:),allocatable :: x
		integer,dimension(1) :: idum
		double precision,dimension(1) :: ddum

		n = size(F,1)
		nnz = size(Kcsr,1)
		nrhs = 1
		maxfct = 1
		mnum = 1
		allocate(x(n))
		allocate(iparm(64))		!set up pardiso control parameter
		do i=1,64
		   iparm(i) = 0
		end do 
		iparm(1) = 1 ! no solver default
		iparm(2) = 2 ! fill-in reordering from metis
		iparm(4) = 0 ! no iterative-direct algorithm
		iparm(5) = 0 ! no user fill-in reducing permutation
		iparm(6) = 0 ! =0 solution on the first n compoments of x
		iparm(8) = 2 ! numbers of iterative refinement steps
		iparm(10) = 13 ! perturbe the pivot elements with 1e-13
		iparm(11) = 1 ! use nonsymmetric permutation and scaling mps
		iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). try iparm(13) = 1 in case of inaccuracy
		iparm(14) = 0 ! output: number of perturbed pivots
		iparm(18) = -1 ! output: number of nonzeros in the factor lu
		iparm(19) = -1 ! output: mflops for lu factorization
		iparm(20) = 0 ! output: numbers of cg iterations

		error  = 0 ! initialize error flag
		msglvl = 0 ! 0=no output, 1=print statistical information
		mtype  = 11 ! real and unsymmetric matrix

		! Initiliaze the internal solver memory pointer. This is only necessary for the first call of the pardiso solver.

		allocate (pt(64))
		do i=1,64
		   pt(i)%dummy =  0 
		end do

		phase = 11 ! Only reordering and symbolic factorization

		call pardiso (pt, maxfct, mnum, mtype, phase, n, Kcsr, ia, ja, &
		idum, nrhs, iparm, msglvl, ddum, ddum, error)
		
!		write(*,*) 'reordering completed ... '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		end if
!		write(*,*) 'number of nonzeros in factors = ',iparm(18)
!		write(*,*) 'number of factorization mflops = ',iparm(19)

		phase = 22 ! only factorization
		call pardiso (pt, maxfct, mnum, mtype, phase, n, Kcsr, ia, ja, &
		idum, nrhs, iparm, msglvl, ddum, ddum, error)
!		write(*,*) 'factorization completed . '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif

		! back substitution and iterative refinement
		iparm(8) = 2 ! max numbers of iterative refinement steps
		phase = 33 ! only solving
		call pardiso (pt, maxfct, mnum, mtype, phase, n, Kcsr, ia, ja, &
		idum, nrhs, iparm, msglvl, F, x, error)
		write(*,*) 'solve completed ... '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif
			  
		1000 continue
		! termination and release of memory
		phase = -1 ! release internal memory
		call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
		idum, nrhs, iparm, msglvl, ddum, ddum, error1)

		if (error1 /= 0) then
		   write(*,*) 'the following error on release stage was detected: ', error1
		   stop 1
		endif

		if ( error /= 0 ) stop 1

	end subroutine solvesystem

end module conductionfem
