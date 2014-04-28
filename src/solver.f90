include 'mkl_pardiso.f90'	! Added for pardiso

module solver

	implicit none
	contains

!-------------------------------------------------------------------
!	Minimal Residual Method
!-------------------------------------------------------------------
	subroutine minres(acsr,ia,ja,b,maxiter,x,iter)
		integer :: maxiter,iter,i,j,k,n,ia(:),ja(:)
		real(8),parameter :: cc=1e-7
		real(8) :: alpha,acsr(:),b(:)
		real(8),allocatable :: x(:),r(:),p(:),temp(:)

		n = size(b,1)
		allocate(x(n))
		allocate(r(n))
		allocate(p(n))
		allocate(temp(n))

		x = 0.d0
		call mkl_dcsrgemv("N",n,acsr,ia,ja,x,p)
		r = b-p
		call mkl_dcsrgemv("N",n,acsr,ia,ja,r,p)

		do i=1,maxiter
			if(mod(i,1000).eq.0) then
				print *, "iteration: ", i
				print *, "residual ratio: ", norm2(r)/cc
			end if
			alpha = dot_product(p,r)/dot_product(p,p)
			x = x + alpha*r
			r = r-alpha*p
			if(norm2(r) .lt. cc) then
				iter = i
				return
			end if
			if(i.eq.maxiter) then
				print *, "Maximum iterations reached."
				print *, "Convergence not achieved."
				print *, "Norm of residual: ", norm2(r)
				print *, "Convergence criterion: ", cc
				if((norm2(r)/cc) .lt. 2.d0) then
					print *, "The residual is within a small ",		&
					"range of the convergence criterion."
					print *, "Perhaps increasing iteration ",		&
					"count may help."
				end if
			end if
			call mkl_dcsrgemv("N",n,acsr,ia,ja,r,p)
		end do

	end subroutine minres

!-------------------------------------------------------------------
!	END Minimal Residual Method
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	Residual Norm Method
!-------------------------------------------------------------------
	subroutine resnorm(acsr,ia,ja,b,maxiter,x,iter)
		integer :: maxiter,iter,i,j,k,n,ia(:),ja(:)
		real(8),parameter :: cc=0.001d0
		real(8) :: alpha,acsr(:),b(:)
		real(8),allocatable :: x(:),r(:),p(:),temp(:)

		n = size(b,1)
		allocate(x(n))
		allocate(r(n))
		allocate(p(n))
		allocate(temp(n))

		x = 0.d0
		call mkl_dcsrgemv("N",n,acsr,ia,ja,x,p)
		r = b-p

		do i=1,maxiter
			if(mod(i,1000).eq.0) then
				print *, "iteration: ", i
				print *, "residual ratio: ", norm2(r)/cc
			end if
			call mkl_dcsrgemv("T",n,acsr,ia,ja,r,p)
			call mkl_dcsrgemv("N",n,acsr,ia,ja,p,temp)
			alpha = (norm2(p)**2.d0)/(norm2(temp)**2.d0)
			x = x + alpha*p
			if(mod(i,10).eq.0) then
				call mkl_dcsrgemv("N",n,acsr,ia,ja,x,temp)
				r = b-temp
			else
				r = r - alpha*temp
			end if
			if(norm2(r) .lt. cc) then
				iter = i
				return
			end if
			if(i.eq.maxiter) then
				print *, "Maximum iterations reached."
				print *, "Convergence not achieved."
				print *, "Norm of residual: ", norm2(r)
				print *, "Convergence criterion: ", cc
				if((norm2(r)/cc) .lt. 2.d0) then
					print *, "The residual is within a small ",		&
					"range of the convergence criterion."
					print *, "Perhaps increasing iteration ",		&
					"count may help."
				end if
			end if
		end do
	end subroutine resnorm
!-------------------------------------------------------------------
!	END Residual Norm Method
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	BiConjugate Gradient (Stabilised) Method
!-------------------------------------------------------------------
	subroutine bicgstab(acsr,ia,ja,b,maxiter,x,iter)
		integer :: maxiter,iter,i,j,k,n,ia(:),ja(:)
		real(8),parameter :: cc=1e-9
		real(8) :: alpha,beta,delta0,delta,delta_old,omega,			&
		acsr(:),b(:)
		real(8),allocatable :: x(:),r(:),p(:),s(:),rst(:),			&
		temp1(:),temp2(:)

		n = size(b,1)
		allocate(x(n))
		allocate(r(n))
		allocate(p(n))
		allocate(s(n))
		allocate(rst(n))
		allocate(temp1(n))
		allocate(temp2(n))

		x = 0.d0
		call mkl_dcsrgemv("N",n,acsr,ia,ja,x,temp1)
		r = b-temp1
		rst = 1.d0
		p = r
		delta = dot_product(rst,r)
		write(*,'(a,1x,f15.3)') "Starting delta: ", delta
		delta0 = delta

		do i=1,maxiter
			if(mod(i,1000).eq.0) then
				write(*,'(a,1x,i6)') 'Iteration number: ',i
				write(*,'(a,1x,f15.3)') "Residual ratio: ", norm2(r)/cc
			end if
			call mkl_dcsrgemv("N",n,acsr,ia,ja,p,temp1)	! temp1=A*p
			alpha = delta/dot_product(rst,temp1)
			s = r - alpha*temp1
			call mkl_dcsrgemv("N",n,acsr,ia,ja,s,temp2)	! temp2=A*s
			omega = dot_product(s,temp2)/dot_product(temp2,temp2)
			x = x + alpha*p + omega*s
			r = s - omega*temp2
			delta_old = delta
			delta = dot_product(rst,r)
			beta = (delta/delta_old)*(alpha/omega)
			p = r + beta*(p - omega*temp1)
			if(norm2(r) .lt. cc) then
				iter = i
				return
			end if
			if(i.eq.maxiter) then
				write(*,'(a)') "Maximum iterations reached."
				write(*,'(a)') "Convergence not achieved."
				write(*,'(a,1x,f15.3)') "Norm of residual: ", norm2(r)
				write(*,'(a,1x,f15.3)') "Convergence criterion: ", cc
				if((norm2(r)/cc) .lt. 2.d0) then
					write(*,'(a)') "The residual is within a small",&
					"range of the convergence criterion."
					write(*,'(a)') "Perhaps increasing iteration ",	&
					"count may help."
				end if
			end if
		end do
	end subroutine bicgstab
!-------------------------------------------------------------------
!	END BiConjugate Gradient (Stabilised) Method
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!	PARDISO Direct Solver
!-------------------------------------------------------------------
	subroutine solvesystem(Kcsr,F,ia,ja,x)
		use mkl_pardiso
		integer,parameter :: dp = kind(1.0d0)
		type(mkl_pardiso_handle),dimension(:),allocatable  :: pt
		integer :: i,maxfct,mnum,mtype,phase,n,nrhs,error,msglvl,	&
		nnz,error1
		integer,dimension(:),allocatable :: iparm
		integer,dimension(:),intent(in) :: ia,ja
		real(8),dimension(:) :: Kcsr,F
		real(8),dimension(:),allocatable :: x
		integer,dimension(1) :: idum
		real(8),dimension(1) :: ddum

		n = size(F,1)
		nnz = size(Kcsr,1)
		nrhs = 1
		maxfct = 1
		mnum = 1

		if(not(allocated(x)))	allocate(x(n))
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
		iparm(13) = 0 ! maximum weighted matching algorithm is
					  !switched-off (default for symmetric).
					  ! try iparm(13) = 1 in case of inaccuracy
		iparm(14) = 0 ! output: number of perturbed pivots
		iparm(18) = -1 ! output: number of nonzeros in the factor lu
		iparm(19) = -1 ! output: mflops for lu factorization
		iparm(20) = 0 ! output: numbers of cg iterations

		error  = 0 ! initialize error flag
		msglvl = 0 ! 0=no output, 1=print statistical information
		mtype  = 11 ! real and unsymmetric matrix

		! Initiliaze the internal solver memory pointer.
		! This is only necessary for the first call of the solver.

		allocate (pt(64))
		do i=1,64
		   pt(i)%dummy =  0
		end do

		phase = 11 ! Only reordering and symbolic factorization

		call pardiso (pt,maxfct,mnum,mtype,phase,n,Kcsr,ia,ja, 		&
		idum, nrhs, iparm, msglvl, ddum, ddum, error)

		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		end if

		phase = 22 ! only factorization
		call pardiso (pt,maxfct,mnum,mtype,phase,n,Kcsr,ia,ja, 		&
		idum, nrhs, iparm, msglvl, ddum, ddum, error)
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif

		! back substitution and iterative refinement
		iparm(8) = 2 ! max numbers of iterative refinement steps
		phase = 33 ! only solving
		call pardiso (pt,maxfct,mnum,mtype,phase,n,Kcsr,ia,ja, 		&
		idum, nrhs, iparm, msglvl, F, x, error)
		write(*,*) 'solve completed ... '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif

		1000 continue
		! termination and release of memory
		phase = -1 ! release internal memory
		call pardiso (pt,maxfct,mnum,mtype,phase,n,ddum,idum,idum, 	&
		idum, nrhs, iparm, msglvl, ddum, ddum, error1)

		if (error1 /= 0) then
		   write(*,*) 'the following release error was detected: ',	&
		   error1
		   stop 1
		endif

		if ( error /= 0 ) stop 1

	end subroutine solvesystem
!-------------------------------------------------------------------
!	END PARDISO Direct Solver
!-------------------------------------------------------------------

end module solver
