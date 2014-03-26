module transient

	use assembly

	implicit none
	contains

	subroutine getcapacitance(elDom,elVol,DomCs,DomRhos,elCp)
		integer :: elDom
		real(8) :: elVol,DomCs(:),DomRhos(:),elCp(4,4)

		elCp = 1.d0;
		do i=1,4
			elCp(i,i) = 2.d0
		end do
		elCp = DomCs(elDom)*DomRhos(elDom)*(elVol/120.d0)*elCp
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
	syInit)
		integer :: i,j,k,n,nstep,iaSt(:),jaSt(:),iaCp(:),jaCp(:)
		real(8) :: tstep,tfinal,sySt(:),syCp(:),sySrc(:),syInit(:)
		real(8),allocatable :: kTprod(:),rhs(:),diffT(:),newT(:),	&
		k1(:),k2(:),k3(:),k4(:)

		n = size(sySrc,1)
		allocate(kTprod(n))
		allocate(rhs(n))
		allocate(diffT(n))
		allocate(newT(n))
		allocate(k1(n))
		allocate(k2(n))
		allocate(k3(n))
		allocate(k4(n))

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

	end subroutine transientsolve

!-------------------------------------------------------------------
!	Subroutine to get time gradient at each initial step of RK
!-------------------------------------------------------------------

	subroutine timegradient(syCp,iaCp,jaCp,rhs,diffT)
		integer :: iter,iaCp(:),jaCp(:)
		real(8) :: syCp(:),rhs(:),diffT(:)

		call bicgstab(syCp,iaCp,jaCp,rhs,100000,diffT,iter)

	end subroutine timegradient

end module transient
