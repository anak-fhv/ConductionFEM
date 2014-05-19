module transient

	use assembly
	use solver
	use postproc

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
	useRK,syInit,noVerts,connTab,nDoms,doElems)
		integer :: nDoms,iaSt(:),jaSt(:),iaCp(:),jaCp(:),doElems(:),&
		connTab(:,:)
		real(8) :: sySt(:),syCp(:),sySrc(:),syInit(:),noverts(:,:)
		logical :: useRK

		if(useRK) then
			call transientRK(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,	&
			syInit,noVerts,connTab,nDoms,doElems)
		else
			call transientFD(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,	&
			syInit,noVerts,connTab,nDoms,doElems)
		end if
	end subroutine transientsolve

	subroutine transientFD(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,		&
	syInit,noVerts,connTab,nDoms,doElems)
		integer,parameter :: trfileno=888
		integer :: i,j,fno,n,nv,ntstep,nDoms,iter,iaSt(:),jaSt(:),	&
		iaCp(:),jaCp(:),doElems(:),connTab(:,:)
		real(8) :: theta,tstep,tfinal,sySt(:),syCp(:),sySrc(:),		&
		syInit(:),noVerts(:,:)
		real(8),allocatable :: CKLhs(:),CKRhs(:),FRhs(:),Tnew(:),	&
		kTRhs(:),initGuess(:)
		character(*),parameter :: objdir="../obj/",fres="res",		&
		fext=".vtk"
		character(len=100) :: fName,sysCall

		n = size(sySrc,1)
		allocate(FRhs(n))
		allocate(kTRhs(n))
		allocate(initGuess(n))
		nv = size(sySt,1)
		allocate(CKLhs(nv))
		allocate(CKRhs(nv))

		tstep = 0.001d0
		tfinal= 2.d0
		ntstep = tfinal/tstep + 1
		theta = 1.d0
		CKLhs = theta*tstep*sySt + syCp
		CKRhs = syCp - (1.d0-theta)*tstep*sySt
		initGuess = 0.d0

		do i=1,ntstep
			FRhs = theta*sySrc + (1.d0-theta)*sySrc
			call mkl_dcsrgemv("N",n,CKRhs,iaSt,jaSt,syInit,kTRhs)
			FRhs = tstep*FRhs + kTRhs
			call bicgstab(CKLhs,iaCp,jaCp,FRhs,100000,initGuess,	&
			Tnew,iter)
			syInit = Tnew
			deallocate(Tnew)
!			open(trfileno,file=objdir//ftr)
!			do j=1,n
!				write(trfileno,'(3(f9.4,2x),f9.4)') noVerts(j,:),	&
!				syInit(j)
!			end do
!			close(trfileno)
			if(mod(i,10).eq.0) then
				fno = (i+1)/10
				call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
				syInit)
				write(fName,*) fno
				write(fName,*) objdir//fres//trim(adjustl(fName))//fext
				write(sysCall,*)"mv "//objdir//fres//fext//" "//trim(adjustl(fName))
				write(*,*) sysCall
				call system(sysCall)
			end if
		end do
		
	end subroutine transientFD

	subroutine transientRK(sySt,iaSt,jaSt,syCp,iaCp,jaCp,sySrc,		&
	syInit,noVerts,connTab,nDoms,doElems)
		integer,parameter :: trfileno=888
		integer :: i,j,k,n,ntstep,nDoms,iaSt(:),jaSt(:),iaCp(:),	&
		jaCp(:),doElems(:),connTab(:,:)
		real(8) :: tstep,tfinal,sySt(:),syCp(:),sySrc(:),syInit(:),	&
		noVerts(:,:)
		real(8),allocatable :: kTprod(:),rhs(:),diffT(:),newT(:),	&
		k1(:),k2(:),k3(:),k4(:)
		character(*),parameter :: objdir="../obj/",fres="res",		&
		fext=".vtk"
		character(len=100) :: fName,sysCall

		n = size(sySrc,1)
		allocate(kTprod(n))
		allocate(rhs(n))
		allocate(newT(n))

		tstep = 0.0001
		tfinal = 2
		ntstep = tfinal/tstep + 1

		do i = 1,ntstep
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
			if(allocated(k1)) deallocate(k1)
			if(allocated(k2)) deallocate(k2)
			if(allocated(k3)) deallocate(k3)
			if(allocated(k4)) deallocate(k4)
			call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
			syInit)
			write(fName,*) i
			write(fName,*) objdir//fres//trim(adjustl(fName))//fext
			write(sysCall,*)"mv "//objdir//fres//fext//" "//trim(adjustl(fName))
			write(*,*) sysCall
			call system(sysCall)

!			open(trfileno,file=objdir//ftr)
!			do j=1,n
!				write(trfileno,'(3(f9.4,2x),f9.4)') noVerts(j,:),	&
!				syInit(j)
!			end do
!			close(trfileno)

		end do

	end subroutine transientRK

!-------------------------------------------------------------------
!	Subroutine to get time gradient at each initial step of RK
!-------------------------------------------------------------------

	subroutine timegradient(syCp,iaCp,jaCp,rhs,diffT)
		integer :: n,iter,iaCp(:),jaCp(:)
		real(8) :: syCp(:),rhs(:)
		real(8),allocatable :: initGuess(:),diffT(:)

		n = size(rhs,1)
		allocate(initGuess(n))
		initGuess = 0.d0

		call bicgstab(syCp,iaCp,jaCp,rhs,100000,initGuess,diffT,iter)

	end subroutine timegradient

end module transient
