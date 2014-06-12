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

	subroutine writeresultsvtk(noVerts,connTab,nDoms,doElems,reVals)
		integer,parameter :: fid = 246
		integer :: i,nNodes,nElems,nDoms,nCorners,tetType,			&
		doElems(:),connTab(:,:)
		real(8) :: reVals(:),noVerts(:,:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"res.vtk"

		nNodes = size(noVerts,1)
		nElems = size(connTab,1)
		nCorners = 4
		tetType = 10
		open(fid,file=resfile)
		write(fid,'(a)')"# vtk DataFile Version 1.0"
		write(fid,'(a)')"3D Unstructured Grid of Linear Tetrahedrons"
		write(fid,'(a)')"ASCII"
		write(fid,'(a)')""
		write(fid,'(a)')"DATASET UNSTRUCTURED_GRID"
		write(fid,'(a,2x,i8,2x,a)')"POINTS ",nNodes," double"
		do i=1,nNodes
			write(fid,*) noVerts(i,:)
		end do
		write(fid,'(a)')""
		write(fid,'(a,2x,i8,2x,i8)')"CELLS ",nElems,5*nElems
		do i=1,nElems
			write(fid,'(i2,2x,3(i8,2x),i8)')nCorners,connTab(i,:)-1
		end do
		write(fid,'(a)')""
		write(fid,'(a,2x,i8)')"CELL_TYPES ",nElems
		do i=1,nElems
			write(fid,'(i4)')tetType
		end do
		if(nDoms.gt.1) then
			write(fid,'(a,2x,i8)')"CELL_DATA",nElems
			write(fid,'(a,2x,i2)')"FIELD FieldData", 1
			write(fid,'(a,2x,i2,2x,i8,2x,a)')"Material",1,nElems,"int"
			write(fid,'(5(i4,2x))')doElems-1
		end if
		write(fid,'(a)')""
		write(fid,'(a,2x,i8)')"POINT_DATA ",nNodes
		write(fid,'(a)')"SCALARS temperature double"
		write(fid,'(a)')"LOOKUP_TABLE default"
		do i=1,nNodes
			write(fid,*) reVals(i)
		end do

		close(fid)

	end subroutine writeresultsvtk

	subroutine checkemissiondifference(noVerts,connTab,sfElems,		&
	emFc,reVals,faceEmDiffs)
		integer :: i,j,k,nElems,nNodes,emFc,bFcNo(3),connTab(:,:),	&
		sfElems(:,:)
		real(8) :: aC,T1,T2,T3,em1,em2,reVals(:),noVerts(:,:),		&
		faceEmDiffs(:,:)

		nElems = size(connTab,1)
		nNodes = size(noVerts,1)

		do i=1,nElems
			if(any(sfElems(i,:)) == emFace) then
				do j=1,4
					if(sfElems(i,j) == emFace) then
						call bfacenodes(j,bFcNo)
						bFcNo = connTab(i,bFcNo)
						fcA = facearea(noVerts(bFcNo,:))
						T1 = reVals(bFcNo(1))
						T2 = reVals(bFcNo(2))
						T3 = reVals(bFcNo(3))
						em1 = fcA*sigB*aC/(30.d0*(T3-T1))
						em1 = em1*((T2**6.d0-T3**6.d0)/(T2-T3) + 	&
						(T2**6.d0-T1**6.d0)/(T1-T2))
						em2 = fcA*sigB*aC*((T1+T2+T3)/3.d0)**4
						faceEmDiffs(i,:) = (/em1,em2,em1-em2/)
					end if
				end do
			end if
		end do
	end subroutine checkemissiondifference

end module postproc
