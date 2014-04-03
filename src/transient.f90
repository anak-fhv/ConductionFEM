module transient

	use assembly
	use solver

	implicit none
	contains

	subroutine getcapacitance(elDom,elVol,domCs,domRhos,elCp)
		integer :: i,elDom
		real(8) :: elVol,domCs(:),domRhos(:),elCp(4,4)

		elCp = 1.d0;
		do i=1,4
			elCp(i,i) = 2.d0
		end do
		! Note: 120 is needed because elVol = 6V has been used
		elCp = domCs(elDom)*domRhos(elDom)*(elVol/120.d0)*elCp
	end subroutine getcapacitance

	subroutine assemblecapacitance(noElemCap,elNodes,elCp)
		type(noderow) :: noElemCap(4)
		integer :: elNodes(4)
		real(8) :: elCp(4,4)

		call assemble_noderows(noElemCap,elNodes,elCp)
	end subroutine assemblecapacitance

!-------------------------------------------------------------------
!	Very important note:
!	When transientsolve is called, the subroutine for setting up
!	the final equations, called by the regular FEM routine,shall 
!	have to be discarded. This is because the stiffness matrix, as
!	well as the source term, need to be preserved in their original
!	forms for the transient simulation.
!-------------------------------------------------------------------

	subroutine transientsolve(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,	&
	useRK,syInit,noVerts)
		integer :: iaSt(:),jaSt(:),iaCp(:),jaCp(:)
		real(8) :: sySt(:),syCp(:),sySrc(:),syInit(:),noverts(:,:)
		logical :: useRK

		useRK = .false.
		if(useRK) then
			call transientRK(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,syInit)
		else
			call transientFD(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,syInit,noVerts)
		end if
	end subroutine transientsolve

	subroutine transientFD(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,		&
	syInit,noVerts)
		integer,parameter :: trfileno=888,indices=(/1,3,5,7/)
		integer :: i,j,fno,n,nv,ntstep,iter,iaSt(:),jaSt(:),iaCp(:),&
		jaCp(:)
		real(8) :: theta,tstep,tfinal,sySt(:),syCp(:),sySrc(:),		&
		syInit(:),noVerts(:,:)
		real(8),allocatable :: CKLhs(:),CKRhs(:),FRhs(:),Tnew(:),	&
		kTRhs(:)
		character(*),parameter :: objdir="../obj/",ftr="trans.out"
		character(6) :: fname

		n = size(sySrc,1)
		allocate(FRhs(n))
		allocate(kTRhs(n))
		nv = size(sySt,1)
		allocate(CKLhs(nv))
		allocate(CKRhs(nv))

		tstep = 0.0001d0
		tfinal= 2.d0
		nstep = tfinal/tstep + 1
		theta = 1.d0
		CKLhs = theta*tstep*sySt + syCp
		CKRhs = syCp - (1.d0-theta)*tstep*sySt

		open(trfileno,file=objdir//ftr)
		do i=1,ntstep
			FRhs = theta*sySrc + (1.d0-theta)*sySrc
			call mkl_dcsrgemv("N",n,CKRhs,iaSt,jaSt,syInit,kTRhs)
			FRhs = tstep*FRhs + kTRhs
			call bicgstab(CKLhs,iaCp,jaCp,FRhs,100000,Tnew,iter)
			write(trfileno,'(3(f9.4,2x),f9.4)') Tnew(indices)
			syInit = Tnew
			deallocate(Tnew)
		end do
		close(trfileno)
	end subroutine transientFD

	subroutine transientRK(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,		&
	syInit)
		integer,parameter :: trfileno=888,indices=(/1,3,5,7/)
		integer :: i,j,k,fno,n,nstep,iaSt(:),jaSt(:),iaCp(:),jaCp(:)
		real(8) :: tstep,tfinal,sySt(:),syCp(:),sySrc(:),syInit(:)
		real(8),allocatable :: kTprod(:),rhs(:),diffT(:),newT(:),	&
		k1(:),k2(:),k3(:),k4(:)
		character(*),parameter :: objdir="../obj/",ftr="trans.out"

		n = size(sySrc,1)
		allocate(kTprod(n))
		allocate(rhs(n))
		allocate(newT(n))

		tstep = 0.0001
		tfinal = 2
		nstep = tfinal/tstep + 1

		open(trfileno,file=objdir//ftr)
		do i = 1,nstep
			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,syInit,kTprod)
			rhs = sySrc - kTprod
			call timegradient(syCp,iaCp,jaCp,rhs,k1)
			newT = syInit + (tstep/2.d0)*k1
			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
			rhs = sySrc - kTprod
			call timegradient(syCp,iaCp,jaCp,rhs,k2)
			newT = syInit + (tstep/2.d0)*k2
			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
			rhs = sySrc - kTprod
			call timegradient(syCp,iaCp,jaCp,rhs,k3)
			newT = syInit + tstep*k3
			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
			rhs = sySrc - kTprod
			call timegradient(syCp,iaCp,jaCp,rhs,k4)
			syInit = syInit + (tstep/6)*(k1 + 2*k2 + 2*k3 + k4)

			write(trfileno,'(3(f9.4,2x),f9.4)') syInit(indices)
			if(allocated(k1)) deallocate(k1)
			if(allocated(k2)) deallocate(k2)
			if(allocated(k3)) deallocate(k3)
			if(allocated(k4)) deallocate(k4)
		end do

	end subroutine transientRK

!-------------------------------------------------------------------
!	Subroutine to get time gradient at each initial step of RK
!-------------------------------------------------------------------

	subroutine timegradient(syCp,iaCp,jaCp,rhs,diffT)
		integer :: iter,iaCp(:),jaCp(:)
		real(8) :: syCp(:),rhs(:)
		real(8),allocatable :: diffT(:)

		call bicgstab(syCp,iaCp,jaCp,rhs,100000,diffT,iter)

	end subroutine timegradient

end module transient
