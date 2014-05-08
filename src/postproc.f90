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

	subroutine writeresultsvtk(noVerts,connTab,revals)
		integer,parameter :: fid = 246
		integer :: nNodes,nElems,nCorners,tetType,connTab(:,:)
		real(8) :: reVals(:),noVerts(:,:)
		character(*),parameter :: objdir = "../obj/",				&
								  resfile = objdir//"res.vtk"

		nNodes = size(noVerts,1)
		nElems = size(connTab,1)
		nCorners = 4
		tetType = 10
		open(fid,file=resfile)
		write(fid,*)"# vtk DataFile Version 1.0"
		write(fid,*)"3D Unstructured Grid of Linear Tetrahedrons"
		write(fid,*)"ASCII"
		write(fid,*)""
		write(fid,*)"DATASET UNSTRUCTURED_GRID"
		write(fid,'(a,2x,i8,2x,a)')"POINTS ",nNodes," double"
		do i=1,nNodes
			write(fid,*) noVerts(i,:)
		end do
		write(fid,*)""
		write(fid,'(a,2x,i8,2x,i8)')"CELLS ",nElems,5*nElems
		do i=1,nElems
			write(fid,'(3(i8,2x),i8)')nCorners,connTab(i)
		end do
		write(fid,*)""
		write(fid,*)"CELL_TYPES ",nElems
		do i=1,nElems
			write(fid,'(i4)')tetType
		end do
		write(fid,*)""
		write(fid,*)"POINT_DATA ",nNodes
		write(fid,*)"SCALARS temperature double"
		write(fid,*)"LOOKUP_TABLE default"
		do i=1,nNodes
			write(fid,*) revals(i)
		end do

	end subroutine writeresultsvtk

end module postproc
