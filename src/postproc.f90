module postproc

	use element
	implicit none

	contains

	subroutine getflowrates(noVerts,connTab,doElems,domKs,sfElems,	&
	reVals,bl,bh,bindex,qdimlow,qdimhigh)
		integer :: i,j,k,ct1,ct2,bloc,bindex,ax,nElems,bl(:),bh(:),	&
		elnodes(4),fcnodes(3),bfcs(4),doElems(:),connTab(:,:),		&
		sfElems(:,:)
		real(8) :: ev6,elK,qdimlow,qdimhigh,adimlow,adimhigh,		&
		fcarea,fntemp(4),fc(3,3),elVerts(4,3),spfns(4,4),domKs(:),	&
		noVerts(:,:),reVals(:)

		qdimlow = 0.d0
		qdimhigh = 0.d0
		ct1 = 0
		adimlow = 0.d0
		ct2 = 0
		adimhigh = 0.d0

		ax = bindex + 1
		nElems = size(connTab,1)

		do i=1,nElems
			bfcs = sfElems(i,:)
			elNodes = connTab(i,:)
			elVerts = noVerts(elNodes,:)
			elK = domKs(doElems(i))
			call shapefunctions(elVerts,ev6,spfns)
			if(all(bfcs==0)) then
				continue
			else
				do j=1,4
					if(bfcs(j) /= 0) then
						bloc = bfcs(j)
						do k=1,size(bl,1)
							if(bloc == bl(k)) then
								ct1 = ct1+1
								call bfacenodes(j,fcnodes)
								fc = elVerts(fcnodes,:)
								fcarea = facearea(fc)
								adimlow = adimlow + fcarea
								fntemp = spfns(:,ax)
								qdimlow = qdimlow + elK*fcarea*		&
								dot_product(fntemp,reVals(elNodes))
							end if
						end do
						do k=1,size(bh,1)
							if(bloc == bh(k)) then
								ct2 = ct2+1
								call bfacenodes(j,fcnodes)
								fc = elVerts(fcnodes,:)
								fcarea = facearea(fc)
								adimhigh = adimhigh + fcarea
								fntemp = spfns(:,ax)
								qdimhigh = qdimhigh + elK*fcarea*	&
								dot_product(fntemp,reVals(elNodes))
							end if
						end do
					end if
				end do
			end if
		end do

		print *, "Faces low: ", ct1, "Faces high: ",ct2
		print *, "Area on face low: ",adimlow
		print *, "Area on face high: ",adimhigh
		qdimlow = qdimlow/adimlow
		qdimhigh = qdimhigh/adimhigh
		
	end subroutine getflowrates

end module postproc
