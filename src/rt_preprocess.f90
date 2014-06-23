! Raytracing pre-process module

module pre_process

    use rt_elements
    use rt_helper
    
    implicit none
    
    contains
    
    ! wrapper for starting pre-processing. Basically checks whether provided input file exists,
    ! and calls appropriate pre-processing routine. In case file does not exist, it checks 
    ! if the respective other file (either *.msh or *.data) exists and calls respective routine.
    subroutine start_preprocessing
    
		logical :: l
		integer :: alloc_status, io_error, i, n, tmp
        integer :: id        
        ! check if file provided by user exists
        ! first check for *.data file (since it is faster)
        inquire(file=dataFolder//trim(data_fname)//".data", exist=l)
        if (l .eqv. .true.) then
			 write(*,*) 
	         write(*,*) "using "//trim(data_fname)//".data as input"
	         write(*,*)
             call ReadData(trim(data_fname)//".data")             
        else ! use *.msh file since *.data does not exist
	        inquire(file=dataFolder//trim(data_fname)//".msh", exist=l)
	        if (l .eqv. .false.) then
		        write(*,*) trim(data_fname)//" does not exist! (neither as .msh or .data)"
	            stop
            end if
            write(*,*) 
	        write(*,*) "using "//trim(data_fname)//".msh as input"
	        write(*,*)
            call read_mesh_data(trim(data_fname)//".msh")
	        call WriteData(trim(data_fname)//".msh")
        end if    
        
        ! read additional data files
        write(*,*)
        write(*,*) "reading additional data files ..."
        write(*,*)
        if (RT_setup .eq. 'tomo') then    
            
            open(unit=88, file=dataFolder//input, status='old', action='read')
            read(88,'(i10)') n
            
            allocate(source(n,3), stat = alloc_status)
            call check_io_error(io_error,"reading input tomo",88)
            
            do i =1,size(source,1)
! 			    read(88,'(F8.0,1x,F2.0,1x,F13.8)') source(i,1), source(i,2), source(i,3)
                read(88,*) source(i,:) 
! 			    write(*,*) source(i,1), source(i,2), source(i,3)
            end do
            close(unit=88)
!             write(*,*) "read source" 

            
            do i = 1,size(emSurfNames)
                id = maxval((/1:size(emSurfNames)/),1,emSurfNames(i) == emSurf%name)
!                 write(*,*) emSurf(id)%name
                tmp = size(emSurf(id)%elemData,1)
	            allocate(emSurf(id)%value(tmp), stat=alloc_status)
	            call check_alloc_error(alloc_status, "ems%value array tomo")
	            allocate(emSurf(id)%rays(tmp,2), stat=alloc_status)
	            call check_alloc_error(alloc_status, "ems%rays array tomo")
	            emSurf(id)%rays=0
                call CreateEmissionSurfTOMO(emSurf(id))
                write(*,*) emSurf(id)%name, emSurf(id)%originalID, emSurf(id)%power, size(emSurf(id)%rays) 
            end do
            write(*,*)
            
!             stop
            
        elseif (RT_setup .eq. 'led') then
	        
            ! blue spectrum
            open(unit=88, file=dataFolder//"LED_spektrum_blue.data", status = 'old', action='read',iostat = io_error)
            call check_io_error(io_error,"reading led blue",88)
            
            read(88,'(i8)') n
            
			allocate(spectrumB(n,2), stat=alloc_status)
		    call check_alloc_error(alloc_status, "spectrum data")
            
		    do i=1,n
			    read(88,'(e14.6, 1x, e14.6)') spectrumB(i,:)
            end do
            close(unit=88)
            
            ! yellow spectrum
            open(unit=88, file=dataFolder//"LED_spektrum_yellow.data", status = 'old', action='read',iostat = io_error)
            call check_io_error(io_error,"reading led yellow",88)
            
            read(88,'(i8)') n
            
			allocate(spectrumY(n,2), stat=alloc_status)
		    call check_alloc_error(alloc_status, "yellow spectrum data")
            
		    do i=1,n
			    read(88,'(e14.6, 1x, e14.6)') spectrumY(i,:)
            end do
            close(unit=88)  
            
            ! populate emission surface data
            tmp = size(emSurf(1)%elemData,1)
	        allocate(emSurf(1)%value(tmp), stat=alloc_status)
	        call check_alloc_error(alloc_status, "ems%value array led")
	        allocate(emSurf(1)%rays(tmp,2), stat=alloc_status)
	        call check_alloc_error(alloc_status, "ems%rays array led")
		    emSurf(id)%rays=0
		    call CreateEmissionSurfLED(emSurf(1))
            
		end if
		 
    end subroutine start_preprocessing
    
    
    ! main routine which perform reading of input data and calls routines for
    ! creating emission surfaces and populates tetraElement-type
    subroutine read_mesh_data(file_name)

        character(len=*), intent(in)         :: file_name 
        integer                              :: io_error, read_error, alloc_status, i, j, k, n, h
        integer                              :: nVertices, nTetra, nHexa, nPyr, nWedges, nDomain, nSurface
        integer, dimension(:,:), allocatable :: surfData
        integer, dimension(:), allocatable   :: elemDomain
        real                                 :: t1, t2
        integer                              :: nElem, remainElem
        character(len=100)                   :: dName, surfName, line
    
        call cpu_time(t1) ! just out of interest measure time of routine
        
        write(*,*)
        write(*,*) "Start Reading Input File"
        write(*,*) "==========================="
        write(*,*)
    
        ! open file
        open(unit=81, file=dataFolder//file_name, status='old', action='read', iostat=io_error)       
        call check_io_error(io_error,"opening file",81)
    
        ! read first 3 lines, only the 3rd line is of interest
        do n = 1,3
            read(81,'(A)', iostat=read_error) line
            call check_io_error(read_error,"reading first 3 lines",81)
        end do
    
        ! read internally to corresponding integer values
        read(line,'(7(1x,i8))') nVertices, nTetra, nHexa, nPyr, nWedges, nDomain, nSurface
        if (any([nHexa,nPyr,nWedges] /= 0)) write(*,*) "found non-Tetra elements. These will not be considered yet!"
         
        ! allocate and read vertices 
        allocate(vertices(nVertices,3), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertex array")
        do n = 1,nVertices
            read(81,'(3e14.6)', iostat=read_error) vertices(n,:) 
            call check_io_error(read_error,"reading vertex data",81)
        end do
        ! just some status messages
        write(*,*) "read ",nVertices," vertices" 
        
        ! allocate tetraData        
        allocate(tetraData(nTetra), stat=alloc_status)
        call check_alloc_error(alloc_status, "tetraData array")
        
        ! read tetra elements and assign vertices
        do n = 1,nTetra
            read(81,'(4(1x,i8))',iostat=read_error) tetraData(n)%vertexIds
            call check_io_error(read_error,"reading tetra elements",81)
        end do
        write(*,*) "read ",nTetra," tetraeder elements" 
        write(*,*)
        
        ! read elements for each domain        
        do i = 1,nDomain
    
            ! read number of elements with domain and name of domain
            read(81,'(i8,1x,a8)', iostat=read_error) nElem, dName
            
            allocate(elemDomain(nElem), stat=alloc_status)
            call check_alloc_error(alloc_status, "elemDomain array")
        
            ! read data in local array  
            do n = 1,nElem/10
               read(81,'(10(1x,i8))',iostat=read_error) elemDomain(10*n-9:10*n)
               call check_io_error(read_error,"reading domain elements",81)
            end do     
            remainElem = nElem - (nElem/10)*10
            if (remainElem > 0) then
                read(81,"("//achar(ichar('0')+remainElem)//"(1x,i8))",iostat=read_error) elemDomain(10*(nElem/10)+1:nElem)
                call check_io_error(read_error,"reading remainder domain elements",81)
            end if
            
            ! assign data to tetraElement-type
            forall(n = 1:nElem) tetraData(elemDomain(n))%domain = i
            
            deallocate(elemDomain, stat=alloc_status)
            call check_alloc_error(alloc_status, "elemDomain, dealloaction")
                    
!             write(*,'(1x,A5,i7,A)') "read ",nElem," elements for "//dName
!             write(*,'(1x,A43,i7,1x,i2)') "array for elements in domain has now shape ", shape(elemDomain) 
!             write(*,*)
       
        end do     
        
        ! read surface information
        allocate(emSurf(size(emSurfNames)), stat=alloc_status)
        call check_alloc_error(alloc_status, "ems array")
        k = 1
        do i = 1,nSurface
    
            ! read number of elements within surface and name of aurface
            read(81,'(i8,1x,a)', iostat=read_error) nElem, surfName
        
            allocate(surfData(nElem,2), stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData")
            
            ! read data in local array  
            do n = 1, nElem/5
                read(81,'(5(2x,i8,1x,i1))',iostat=read_error) (surfData(j,1), surfData(j,2),j=5*n-4,5*n)
                call check_io_error(read_error,"reading surface data",81)
            end do         
            remainElem = nElem - (nElem/5)*5
            if (remainElem > 0) then
                read(81,"("//achar(ichar('0')+remainElem)//"(2x,i8,1x,i1))",iostat=read_error) (surfData(j,1), surfData(j,2),j=(nElem/5)*5+1,nElem)
            end if 
       
            ! assign data to tetraElement-type
            ! if a face is on the surface, the tetraeder elements reference itself and the
            ! neighboring face is a negative number with the index of the surface
            ! the interface within a domain is ignored
            if (count(surfName == ignoredSurfaces) == 0) then
	            do n = 1, nElem
	                tetraData(surfData(n,1))%neighbors(surfData(n,2),1) = surfData(n,1)
	                tetraData(surfData(n,1))%neighbors(surfData(n,2),2) = -i
	            end do
            end if
            
            ! if the surface is a emission surface
            if (any(surfName == emSurfNames)) then
             
		        ! remember id of emission surface
                emSurf(k)%originalID = i
                
                ! allocate data
                allocate(emSurf(k)%elemData(nElem,2), stat=alloc_status)
                call check_alloc_error(alloc_status, "ems(k)%elemData array")
                
                ! create emission surface
                emSurf(k)%name = surfName
                emSurf(k)%elemData = surfData
                                    
                ! increase counter
                k = k + 1
                
            end if
            
            ! deallocate surfData
            deallocate(surfData, stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData, dealloaction") 
                
        end do
        
        ! close file
        close(unit=81)
        
        ! assign neighbors
        write(*,*)
        write(*,*) "Creating Element Connections"
        write(*,*) "============================="
        write(*,*) 
        call assign_neighbors
        
        ! check for double entries (could be commented)
        open(unit=99, file=dataFolder//"doubledata.res", status='replace', action='write')     
        do i = 1,nTetra
            do n = 1,4
	            j = 0
		        
		        ! count how often the current neighbor is referenced as neighbor in the element
		        ! could be more than once if more than face are on a surface
		        j = j + count((tetraData(i)%neighbors(n,1) == tetraData(i)%neighbors(:,1)) .and. (tetraData(i)%neighbors(n,1) /= i)) - 1
	         
                ! check if zero entry exist
                if (tetraData(i)%neighbors(n,1) == 0) j = j + 1
                	            
                if (j > 0) then
                    write(99,'(9(i8))') i, tetraData(i)%neighbors(:,1), tetraData(i)%neighbors(:,2)      
                    write(*,*) "double entry"
                    write(*,*) tetraData(i)%neighbors(:,1)
                    write(*,*) tetraData(i)%neighbors(:,2)
                    write(*,*)
                    exit
                end if
            end do
        end do
        close(unit = 99)
        
        ! report timing
        call cpu_time(t2)
        write(*,*)
        write(*,*) "Finished Pre-Processing"
        write(*,*) "=============================" 
        write(*,*) "pre-processing took (in Sec.): ", t2-t1
    
    end subroutine read_mesh_data


    ! main routine for assigning neighbors
    subroutine assign_neighbors
       
        integer                            :: alloc_status, n, k, r
        integer, dimension(:), allocatable :: list, id_start, id_end
        integer, dimension(1)              :: idFirstZero, maxElem
        
        ! create initial list
        allocate(list(size(tetraData)),stat=alloc_status)
        call check_alloc_error(alloc_status, "index list")
        maxElem(1) = size(tetraData)
        do n = 1,maxElem(1)
            list(n) = n
        end do
        
        ! create partitions
        allocate(id_start(npart),stat=alloc_status)
        call check_alloc_error(alloc_status, "id_start list")
        allocate(id_end(npart),stat=alloc_status)
        call check_alloc_error(alloc_status, "id_end list")
        r = maxElem(1)/npart
        do n = 1,npart
            id_start(n) = 1 +  (n-1)*r
            id_end(n) = n*r
            if (n == npart) id_end(n) = maxElem(1)
!             write(*,*) n
!             write(*,*) id_start(n)
!             write(*,*) id_end(n)
        end do
        
        write(*,*) "number of partitions ", npart
        write(*,*) "number of valid elements ", maxElem
        
        ! loop over all partitions and perform searches
        do n = 1,npart
!         do n = 1,1
             
            write(*,*)
            write(*,*) "parsing PARTITION number", n
            write(*,*) "searching in list of size", size(list(id_start(n):id_end(n)))  
            ! search within same partition
            call find_neighbors_single_list(list(id_start(n):id_end(n)))            
            
            ! filter list of partition n by putting zeros at the end of the list
            call filter_list(list(id_start(n):id_end(n)),id_start(n), id_end(n))
            if (id_end(n) == 0) exit
            write(*,*) "size of list after search", size(list(id_start(n):id_end(n)))
            
            ! search in all partitions with larger index number
            do k = npart, n+1, -1
            
                call find_neighbors(list(id_start(n):id_end(n)), list(id_start(k):id_end(k)))
                
                ! filter list of partition n
                call filter_list(list(id_start(n):id_end(n)),id_start(n), id_end(n))
                write(*,*) "size of list after search",size(list(id_start(n):id_end(n)))," in list",k
                if (id_end(n) == 0) exit
                
                ! filter list of partition k  
                call filter_list(list(id_start(k):id_end(k)),id_start(k), id_end(k))                
                
            end do
            
        end do
        
        ! checking that all elements are assigned, both numbers should be the same
        if (count(list(:) == 0) /= maxElem(1)) then
            write(*,'(a,i8)') "numbers of elements in index list ", count(list(:) == 0)
            write(*,'(a,i8)') "number of tetra-elements          ", maxElem(1)
            write(*,*)
        end if
        
        ! deallocate lists
        deallocate(list,stat=alloc_status)
        call check_alloc_error(alloc_status, "index list, dealloaction")
        deallocate(id_start,stat=alloc_status)
        call check_alloc_error(alloc_status, "id_start list, dealloaction")
        deallocate(id_end,stat=alloc_status)
        call check_alloc_error(alloc_status, "id_end list, dealloaction")
    
    end subroutine assign_neighbors
    
    
    ! search for neighbors in 2 lists
    ! search routine:
    ! for each element in the 1st list a neighbor is looked for in the 2nd list 
    subroutine find_neighbors(list1, list2)

        integer, intent(inout)            :: list1(:), list2(:)
        integer                           :: i, j, k, counter, face_i, face_j
        logical, dimension(4)             :: mask

        ! loop over elements in list1 and identify neighbors
        do i = 1,size(list1)
            
            ! count how many neighbors element i already has
            counter = count(tetraData(list1(i))%neighbors(:,1) > 0)
                                   
            ! search all elements in list2
            do j= 1,size(list2)
            
                ! if element j has already 4 neighbors cycle
                if (list2(j) == 0) cycle    
                
                ! look which vertices of element i are in element j   
                mask = .false.
                forall(k=1:4) mask(k) = any(tetraData(list1(i))%vertexIds(k) == tetraData(list2(j))%vertexIds(:))
                             
                ! 3 vertices are equal means the elements have a common face
                if (count(mask) == 3) then
                    
                    counter = counter + 1 ! increase number of neighbors for element i               
                    face_i = return_facenumber(mask)  ! determine face which is shared by both elements
                    call add_neighbor(list1(i), list2(j), face_i)                     
                    
                    ! if element j has four neighbors, add zero into list
                    if (minval(tetraData(list2(j))%neighbors(:,1)) > 0) then
                        list2(j) = 0
                    end if
                    
                end if
                
                ! if tetra element i has 4 neighbors break inner do-loop
                if (counter == 4) then
                    list1(i) = 0  ! assign zero to element so it will not be used in the following searches
                    exit
                end if
                
            end do            
            
        end do
        
    end subroutine find_neighbors

    
    ! search for neighbors in single lists
    ! search routine:
    ! for each element in the list a neighbor is looked for in all following elements    
	subroutine find_neighbors_single_list(list)

		integer, intent(inout) :: list(:)
        integer                :: i, j, k, counter, face_i, face_j, n
        logical, dimension(4)  :: mask

        n = size(list) ! size of list        
        ! loop over elements in list and identify neighbors
        do i = 1,n-1
                                                 
            ! if element has already 4 neighbors
            if (list(i) == 0) cycle
            
            ! count how many neighbors element i already has
            counter = count(tetraData(list(i))%neighbors(:,1) > 0)
                                   
            ! search all elements below
            do j= i+1,n
            
                ! if element j has already 4 neighbors
                if (list(j) == 0) cycle
                
                ! look which vertices of element i are in element j  
                mask = .false.
                forall(k=1:4) mask(k) = any(tetraData(list(i))%vertexIds(k) == tetraData(list(j))%vertexIds(:))
                             
                ! 3 vertices are equal, elements are neighbors
                if (count(mask) == 3) then
                    
                    counter = counter + 1   ! increase number of neighbors in element i
                    face_i = return_facenumber(mask) ! get face of element i which is shared by both elements
                    call add_neighbor(list(i), list(j), face_i)
                    
                    ! if element j has four neighbors, add zero into list
                    if (minval(tetraData(list(j))%neighbors(:,1)) > 0) then
                        list(j) = 0
                    end if
                    
                end if
                
                ! if tetra element i has 4 neighbors break inner do-loop
                if (counter == 4) then                
                    list(i) = 0  ! assign zero to element so it will not be used in the following searches
                    exit
                end if
                
            end do            
            
        end do
        
    end subroutine find_neighbors_single_list

    
    ! filter list for zeros and update right index boundary
    subroutine filter_list(list, id_start, id_end)
    
        integer, intent(inout) :: list(:)
        integer, intent(in)    :: id_start
        integer, intent(out)   :: id_end
        integer, dimension(1)  :: idFirstZero
        
        if (count(list == 0) > 0) then
            call putZerosAway(list)
            idFirstZero = minloc(list)
            if (idFirstZero(1) == 1) then
                id_end = 0
            else
                id_end = id_start + idFirstZero(1) - 2
            end if
        end if
    
    end subroutine filter_list
    
    
    ! add neighbor data into tetraElement-type
    subroutine add_neighbor(tetra1, tetra2, face1)
    
        integer, intent(in)   :: tetra1, tetra2, face1
        logical, dimension(4) :: mask
        integer               :: k, face2
        
        mask = .false.
        forall(k=1:4) mask(k) = any(tetraData(tetra2)%vertexIds(k) == tetraData(tetra1)%vertexIds(:))
        face2 = return_facenumber(mask)
                    
        ! assign data in tetraElement-type
        tetraData(tetra1)%neighbors(face1,1) = tetra2
        tetraData(tetra1)%neighbors(face1,2) = face2
        tetraData(tetra2)%neighbors(face2,1) = tetra1
        tetraData(tetra2)%neighbors(face2,2) = face1
    
    end subroutine add_neighbor
    
      
    ! create emission surface from surface data which includes
    ! calculating the area of each face and storing it in a
    ! cumulative manner   
    subroutine CreateEmissionSurfLED(ems)

        type(emissionSurface), intent(inout) :: ems        
        integer                              :: i, alloc_status
        integer, dimension(3)                :: vertIds
        real(dp), dimension(3)               :: p1, p2, p3 
        
        ! it's a single emission surface with constant temperature
	    do i =1,size(ems%value)
                
           ! determine vertices for calculating the face area
           call return_facevertIds(ems%elemData(i,2),vertIds)   
           call return_coords(tetraData(ems%elemData(i,1)), vertIds, p1, p2, p3)                        
           if (i == 1) then
                ems%value(1) = 0.5_dp*norm(cross(p2-p1,p3-p1))
           else 
                ems%value(i) = 0.5_dp*norm(cross(p2-p1,p3-p1)) + ems%value(i-1)
           end if 
                
	    end do
        
	    ! normalize value
	    ems%totalvalue = ems%value(size(ems%value))
	    ems%value = ems%value/ems%totalvalue
	    ems%power = Eemitted
        
    end subroutine CreateEmissionSurfLED
    
    subroutine CreateEmissionSurfTOMO(ems)

        type(emissionSurface), intent(inout) :: ems        
        integer                              :: i, alloc_status
        integer, dimension(3)                :: vertIds
        real(dp), dimension(3)               :: p1, p2, p3 
        real(dp)                             :: temp 
        
        ! check kind of emission surface
        if (any(ems%name == ignoredSurfaces) .eqv. .false.) then
	        
	        ! get temperature at the wall
	        if (ems%name .eq. emSurfNames(1)) then
		        temp = tbounds(1)
		    elseif (ems%name .eq. emSurfNames(2)) then
			    temp = tbounds(2)
	        else
		        write(*,*) "something went wrong in chossing wall temperature!"
		        write(*,*) ems%name
		        write(*,*) emSurfNames
		        stop
	        end if  
	        
		    ! it is one of the wall boundaries
	        do i = 1,size(ems%value)    
			    ! determine vertices for calculating the face area
		        call return_facevertIds(ems%elemData(i,2),vertIds)   
	            call return_coords(tetraData(ems%elemData(i,1)), vertIds, p1, p2, p3)                        
	            if (i == 1) then
		            ems%value(1) = 0.5_dp*norm(cross(p2-p1,p3-p1))*sbk*temp**4
	            else 
	                ems%value(i) = 0.5_dp*norm(cross(p2-p1,p3-p1))*sbk*temp**4 + ems%value(i-1)
	            end if
		    end do
		    
        else
	        
	        ! it is the interface surface
            do i = 1,size(ems%value)
                if (i == 1) then
		            ems%value(1) = alpha*source(i,3)
		        else
			        ems%value(i) = alpha*source(i,3) + ems%value(i-1)
			    end if
			    
			    ! check if anything is good
			    if ((ems%elemData(i,1) .ne. source(i,1)) .or. (ems%elemData(i,2) .ne. source(i,2))) then
				    write(*,*) source(i,:)
				    write(*,*) ems%elemData(i,:)
				    stop
				end if   
			end do
        end if
        
	    ! normalize area
	    ems%power = ems%value(size(ems%value))
	    ems%value = ems%value/ems%power
        
    end subroutine CreateEmissionSurfTOMO
    
    
    ! write out data from pre-processing
    subroutine WriteData(file_name)
    
        character(len=*), intent(in)      :: file_name
        character(len=100)                :: fname
        integer                           :: io_error, n, k, write_error
        
        ! create file from exisiting file-name
        fname = dataFolder//file_name(1:index(file_name,".msh")-1)//".data"
        open(unit=82, file=fname, status='replace', action='write', iostat=io_error)       
        call check_io_error(io_error,"opening file for write-out",82)
        
        ! write values for number of various elements
        write(82,'(3(1x,i9))') size(vertices,1), size(tetraData), size(emSurf)
        
        ! write vertices
        do n = 1,size(vertices,1)
            write(82,'(3e14.6)', iostat=write_error) vertices(n,:) 
            call check_io_error(write_error,"writing vertex data",82)
        end do
        
        ! write tetraData (vertices and domain)
        do n = 1,size(tetraData)
            write(82,'(5(1x,i8))',iostat=write_error) tetraData(n)%vertexIds, tetraData(n)%domain 
            call check_io_error(write_error,"writing tetra data",82)
        end do
        
        ! write tetraData (connections)
        do n = 1, size(tetraData)
            write(82,'(1x,4(i8,1x,i4))',iostat=write_error) transpose(tetraData(n)%neighbors) 
            call check_io_error(write_error,"writing tetra data part 2",82)
        end do
        
        ! write emission surface information
        do n = 1, size(emSurf)
            write(82,'(1x,i8,1x,a)',iostat=write_error) size(emSurf(n)%elemData,1), trim(emSurf(n)%name)
            call check_io_error(write_error,"writing emission surface info",82)
            
            do k=1,size(emSurf(n)%elemData,1)
!                 write(82,'(1x,i8,1x,i2,2x,e14.6)',iostat=write_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2), emSurf(n)%area(k)
                write(82,'(1x,i8,1x,i2)',iostat=write_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2)
                call check_io_error(write_error,"writing emission surface data",82)
            end do
!             write(82,'(e14.6)') emSurf(n)%totalarea
            write(82,'(i4)') emSurf(n)%originalID
        end do    
        
        ! close file
        close(unit=82)
  
    end subroutine WriteData
   
   
    ! read data from exisiting pre-processing file
    subroutine ReadData(file_name)
    
        character(len=*), intent(in) :: file_name
        integer                      :: io_error, n, k, read_error, nVertices, nTetra, nSurf, alloc_status, nElem                                        
        
        ! open file
        open(unit=83, file=dataFolder//file_name, status='old', action='read', iostat=io_error)       
        call check_io_error(io_error,"opening file for reading pre-processed data",83)
        
        ! read information about element sizes
        read(83,'(3(1x,i9))') nVertices, nTetra, nSurf
        
        ! allocate vertices and read vertices
        allocate(vertices(nVertices,3), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertex array")
        do n = 1,nVertices
            read(83,'(3e14.6)', iostat=read_error) vertices(n,:) 
            call check_io_error(read_error,"reading vertex data",83)
        end do
        
        ! allocate tetra data and read information
        allocate(tetraData(nTetra), stat=alloc_status)
        call check_alloc_error(alloc_status, "tetraData array")       
        do n = 1,nTetra
            read(83,'(5(1x,i8))',iostat=read_error) tetraData(n)%vertexIds, tetraData(n)%domain 
            call check_io_error(read_error,"reading tetra data",83)
        end do
        do n = 1, nTetra
            read(83,'(1x,4(i8,1x,i4))',iostat=read_error) (tetraData(n)%neighbors(k,1), tetraData(n)%neighbors(k,2), k=1,4)
            call check_io_error(read_error,"reading tetra data part 2",83)
        end do
        
        ! allocate emission surface and read information
        allocate(emSurf(nSurf), stat=alloc_status)
        call check_alloc_error(alloc_status, "emSurf array")
        do n = 1, nSurf
            read(83,'(1x,i8,1x,a)',iostat=read_error) nElem, emSurf(n)%name
            call check_io_error(read_error,"reading emission surface info",83)
            
!             allocate(emSurf(n)%area(nElem), stat=alloc_status)
!             call check_alloc_error(alloc_status, "emSurf%area array")
            allocate(emSurf(n)%elemData(nElem,2), stat=alloc_status)
            call check_alloc_error(alloc_status, "emSurf%elemData array")
            
            do k=1,nElem
!                 read(83,'(1x,i8,1x,i2,2x,e14.6)',iostat=read_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2), emSurf(n)%area(k)
                read(83,'(1x,i8,1x,i2)',iostat=read_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2)
                call check_io_error(read_error,"reading emission surface data",83)
            end do
!             read(83,'(e14.6)') emSurf(n)%totalarea
            read(83,'(i4)') emSurf(n)%originalID
        end do    
        
        ! close file
        close(unit=83)
  
   end subroutine ReadData

end module pre_process