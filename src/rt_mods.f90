! RayTracing modules
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 26.02.2014
! version: 0.2

module rt_constants

    implicit none
    integer, parameter :: dp=selected_real_kind(p=14)
    
end module rt_constants


module rt_types

    use rt_constants
    implicit none
    
    type :: tetraElement
        real(dp), dimension(3,4) :: vertices    ! coordinates of vertices
        integer, dimension(4,2)  :: neighbors   ! element number of neighrboring elements 
        !real(dp), dimension(4,4) :: shape_funcs ! shape functions
        integer                  :: domain      ! to which domain the tetra belongs
        ! whether its on a surface, which face
        ! domain number
    end type tetraElement
    ! Comments:
    ! * instead of the neighbor one could also specifiy which face the tetra-elements share 
end module rt_types    


module rt_funcs

    use rt_types
    use rt_constants
    implicit none
    
    contains
        
    integer function return_facenumber(mask)
            
        implicit none
        logical, dimension(4) :: mask
        
        if (all(mask == [.true., .true., .true., .false.])) then
            return_facenumber = 1
        else if (all(mask == [.true., .true., .false., .true.])) then
            return_facenumber = 2
        else if (all(mask == [.false., .true., .true., .true.])) then
            return_facenumber = 3
        else if(all(mask == [.true., .false., .true., .true.])) then
            return_facenumber = 4
        else
            write(*,*) "unknown mask. Can not determine correct face number!"
            stop
        end if
            
    end function return_facenumber
       
end module rt_funcs
            
                
module math_funs
    
    use rt_constants
    implicit none
    
    contains
    
    ! put all zeros at the end of the list
    subroutine putZerosAway(list)
        implicit none
        integer, intent(inout) :: list(:)
        integer                :: nmax, nzeros, value
        integer, dimension(1)  :: id
        
        nmax = size(list)
        nzeros = count(list(:) == 0)
        if (nzeros == 0) return

        do
            id = minloc(list(:)) ! first index of smallest element 
            if (nmax-nzeros == id(1)-1) exit
            value = list(id(1))
            
            ! move zero at the end, and shift remaining list 1 element to the front
            list(id(1):nmax-1) = list(id(1)+1:nmax) 
            list(nmax) =  value
            
        end do
        
    end subroutine putZerosAway

    ! cross product in 3d
    function cross(v1,v2)
        implicit none
        real(dp), dimension(3) :: cross
        real(dp), dimension(3), intent(in)  :: v1,v2
 
        cross(1) = v1(2)*v2(3) - v1(3)*v2(2) 
        cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
 
    end function cross
  
    ! norm 
    function norm(v)
        implicit none
        real(dp)             :: norm
        real(dp), intent(in) :: v(:)
  
        norm = sqrt(dot_product(v,v))
   
   end function norm
  
end module math_funs


module pre_process_data

    use rt_constants
    use rt_types
    use rt_funcs
    use math_funs
    
    implicit none
    
    contains
    
    subroutine read_mesh_data(file_name, npart, tetraData)
    
        implicit none
        character(len=*), intent(in)                               :: file_name 
        integer, intent(in)                                        :: npart
        logical                                                    :: l
        integer                                                    :: io_error, read_error, n, alloc_status, i, j
        character(len=100)                                         :: line  
        integer                                                    :: nVertices, nTetra, nHexa, nPyr, nWedges, nDomain, nSurface
        real(dp), dimension(:,:), allocatable                      :: vertices 
        integer, dimension(:,:), allocatable                       :: surfData
        integer, dimension(:), allocatable                         :: elemDomain
        real                                                       :: t1, t2
        integer                                                    :: nElem, remainElem
        character(len=20)                                          :: dName, surfName
        character(len=1)                                           :: test
        type(tetraElement), dimension(:), allocatable, intent(out) :: tetraData
        integer, dimension(:,:), allocatable                       :: vertIDs
        logical,dimension(:,:), allocatable                        :: flag
    
        call cpu_time(t1) ! just out of interest measure time of routine
        
        ! check if file exists
        inquire(file=file_name, exist=l)
        if (l .eqv. .false.) then
            write (*,*) file_name//" does not exist!"
            stop
        end if
    
        ! open file
        open(unit=21, file=file_name, status='old', action='read', iostat=io_error)       
        call check_io_error(io_error,"opening file")
    
        ! read first 3 lines, only the 3rd line is of interest
        do n = 1,3
            read(21,'(A)', iostat=read_error) line
            call check_io_error(read_error,"reading first 3 lines")
        end do
    
        ! read internally to corresponding integer values
        read(line,'(7(1x,i8))') nVertices, nTetra, nHexa, nPyr, nWedges, nDomain, nSurface
        if (any([nHexa,nPyr,nWedges] /= 0)) write(*,*) "found non-Tetra elements. These will not be considered yet!"
        
    
        ! allocate and read vertices 
        allocate(vertices(3,nVertices), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertex array")
        do n = 1,nVertices
            read(21,'(3(1x,3e14.6))', iostat=read_error) vertices(:,n) 
            call check_io_error(read_error,"reading vertex data")
        end do
        ! just some status messages
        write(*,'(a5,i7,a31,1x,i1,1x,i7)') "read ",nVertices," vertices into array of shape", shape(vertices) 
        write(*,*)
        
        ! allocate tetraData
        allocate(tetraData(nTetra), stat=alloc_status)
        call check_alloc_error(alloc_status, "tetraData array")
        allocate(vertIds(4,nTetra), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertIDs array")
        
        ! read tetra elements and assign vertices
        do n = 1,nTetra
            read(21,'(4(1x,i8))',iostat=read_error) vertIDs(:,n)
            call check_io_error(read_error,"reading tetra elements")
            forall(i=1:4) tetraData(n)%vertices(:,i) = vertices(:,vertIDs(i,n))
        end do
        write(*,'(a5,i7,a)') "read ",nTetra," tetraeder elements" 
        write(*,*)
        
        ! read elements for each domain        
        do i = 1,nDomain
    
            ! read number of elements with domain and name of domain
            read(21,'(i8,1x,a8)', iostat=read_error) nElem, dName
            
            write(*,*) nElem, dName
            
            allocate(elemDomain(nElem), stat=alloc_status)
            call check_alloc_error(alloc_status, "elemDomain array")
        
            do n = 1,nElem/10
               read(21,'(10(1x,i8))',iostat=read_error) elemDomain(10*n-9:10*n)
               call check_io_error(read_error,"reading domain elements")
            end do
        
            remainElem = nElem - (nElem/10)*10
            if (remainElem > 0) then
                read(21,"("//achar(ichar('0')+remainElem)//"(1x,i8))",iostat=read_error) elemDomain(10*(nElem/10)+1:nElem)
                call check_io_error(read_error,"reading remainder domain elements")
            end if
            
            forall(n = 1:nElem) tetraData(elemDomain(n))%domain = i
            
            deallocate(elemDomain, stat=alloc_status)
            call check_alloc_error(alloc_status, "elemDomain, dealloaction")
                    
!             write(*,'(1x,A5,i7,A)') "read ",nElem," elements for "//dName
!             write(*,'(1x,A43,i7,1x,i2)') "array for elements in domain has now shape ", shape(elemDomain) 
!             write(*,*)
       
        end do
        
        ! read surface information
        do i = 1,nSurface
    
            read(21,'(i8,1x,a)', iostat=read_error) nElem, surfName
!             write(*,*) nElem, surfName
        
            allocate(surfData(nElem,2), stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData")
      
            do n = 1, nElem/5
                read(21,'(5(2x,i8,1x,i1))',iostat=read_error) (surfData(j,1), surfData(j,2),j=5*n-4,5*n)
                call check_io_error(read_error,"reading surface data")
            end do
            
            remainElem = nElem - (nElem/5)*5
            if (remainElem > 0) then
                read(21,"("//achar(ichar('0')+remainElem)//"(2x,i8,1x,i1))",iostat=read_error) (surfData(j,1), surfData(j,2),j=(nElem/5)*5+1,nElem)
            end if 
       
            do n = 1, nElem
                tetraData(surfData(n,1))%neighbors(surfData(n,2),1) = i+nTetra
                tetraData(surfData(n,1))%neighbors(surfData(n,2),2) = -1
            end do
            
            deallocate(surfData, stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData, dealloaction") 
                
        end do
     
        
        ! assign neighbors
        call assign_neighbors(vertIDs, npart, tetraData)
        
!         call cpu_time(t2)
!         write(*,'(A,3es14.6)') "post-processing took (in Sec.): ", t2-t1
        
        ! check for double entries
        do i = 1,nTetra
            do n = 1,4
                if ( (tetraData(i)%neighbors(n,1) > nTetra) .or. (tetraData(i)%neighbors(n,1) == 0))  cycle
                j= count(tetraData(i)%neighbors(n,1) == tetraData(i)%neighbors(:,1))  
                if (j > 1) then
                    write(*,*) "double entry"
                    write(*,*) tetraData(i)%neighbors(:,1)
                    write(*,*) tetraData(i)%neighbors(:,2)
                    write(*,*)
                    exit
                end if
            end do
        end do
               
                
        ! dealloaction of some array   
        deallocate(vertices, stat=alloc_status)
        call check_alloc_error(alloc_status, "vertices, dealloaction")
        deallocate(vertIDs, stat=alloc_status)
        call check_alloc_error(alloc_status, "vertIDs, dealloaction")
        
        ! close file
        close(unit=21)
        
        ! report timing
        call cpu_time(t2)
        write(*,'(A,3es14.6)') "post-processing took (in Sec.): ", t2-t1
    
    end subroutine read_mesh_data

    subroutine check_io_error(stat, message)
    
        implicit none
        integer, intent(in)                :: stat
        character(*), intent(in), optional :: message
        
        if ((stat > 0) .and. present(message)) then
            write(*,'(A)') "io-error: "//message
            close(unit=21)
            stop
        else if (stat > 0) then
            write(*,'(A)') "error in input-output/reading/writing procedure!"
            close(unit=21)
            stop
        end if
        
        return
    end subroutine check_io_error
    
    subroutine check_alloc_error(stat, message)
    
        implicit none
        integer, intent(in)                :: stat
        character(*), intent(in), optional :: message
        
        if ((stat /= 0) .and. present(message)) then
            write(*,'(A)') "allocation-error: "//message
            close(unit=21)
            stop
        else if (stat > 0) then
            write(*,'(A)') "error in allocation procedure!"
            close(unit=21)
            stop
        end if
        
        return
    end subroutine check_alloc_error
    
    
    subroutine assign_neighbors(vertIDs, npart, tetraData)
    
        implicit none    
        integer, intent(in)                            :: vertIDs(:,:)
        integer, intent(in)                            :: npart
        type(tetraElement), intent(inout)              :: tetraData(:)
        integer                                        :: alloc_status, n, k
        integer, dimension(:), allocatable             :: list, id_start, id_end
        integer, dimension(1)                          :: idFirstZero, maxElem
        
        ! just some checking
!         if (all([allocated(vertIDs), allocated(tetraData)]) .eqv. .false.) then
!             write(*,*) "At least one of the arrays vertIDs or tetraData is not allocated!"
!             stop
!         end if
        
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
        do n = 1,npart
            id_start(n) = 1 +  (n-1)*(maxElem(1)/npart)
            id_end(n) = n*maxElem(1)/npart
            if (n == npart) id_end(n) = maxElem(1)
        end do
        
        write(*,*) "number of partitions ", npart
        write(*,*) "number of valid elements ", maxElem
        write(*,*)      
        
        ! loop over all partitions and perform searches
        do n = 1,npart
        
            ! search within same partition
            call find_neighbors_single_list(vertIDS, list(id_start(n):id_end(n)), tetraData)
            write(*,*) "reduction after search in single list", count(list(id_start(n):id_end(n)) == 0)/real(id_end(n) - id_start(n) + 1)
            write(*,*)
            ! filter list of partition n
            if (count(list(id_start(n):id_end(n)) == 0) > 0) then
                call putZerosAway(list(id_start(n):id_end(n)))
                idFirstZero = minloc(list(id_start(n):id_end(n)))
                if (idFirstZero(1) == 1) then
                    id_end(n) = 0
                    exit ! list is empty
                else
                    id_end(n) = id_start(n) + idFirstZero(1) - 2
                end if
            end if
            
            ! search in all partitions with larger index number
            do k = n+1,npart
            
                call find_neighbors(vertIDS, list(id_start(n):id_end(n)), list(id_start(k):id_end(k)), tetraData)
                write(*,*) "reduction in high level list", count(list(id_start(n):id_end(n)) == 0)/real(id_end(n) - id_start(n) + 1)
                write(*,*) "reduction in low level list ", count(list(id_start(k):id_end(k)) == 0)/real(id_end(k) - id_start(k) + 1)
                write(*,*)
                
                ! filter list of partition n
                if (count(list(id_start(n):id_end(n)) == 0) > 0) then
                    call putZerosAway(list(id_start(n):id_end(n)))
                    idFirstZero = minloc(list(id_start(n):id_end(n)))
                    if (idFirstZero(1) == 1) then
                        id_end(n) = 0
                        exit ! list is empty
                    else
                        id_end(n) = id_start(n) + idFirstZero(1) - 2
                    end if
                end if
                
                ! filter list of partition k                
                if ( (count(list(id_start(k):id_end(k)) == 0) > 0) ) then
                    call putZerosAway(list(id_start(k):id_end(k)))
                    idFirstZero = minloc(list(id_start(k):id_end(k)))
                    id_end(k) = id_start(k) +  idFirstZero(1) - 2
                end if                   
                
            end do
            
        end do
        
        write(*,'(a,i8)') "numbers of elements in index list ",count(list(:) == 0)
        write(*,'(a,i8)') "number of tetra-elements          ", maxElem(1)
        write(*,*)
        
        ! deallocate lists
        deallocate(list,stat=alloc_status)
        call check_alloc_error(alloc_status, "index list, dealloaction")
        deallocate(id_start,stat=alloc_status)
        call check_alloc_error(alloc_status, "id_start list, dealloaction")
        deallocate(id_end,stat=alloc_status)
        call check_alloc_error(alloc_status, "id_end list, dealloaction")
        
!         write(*,*) offset
        !show allocation of neighbors     
!         do n=120,140
!             write(*,'(5(1x,i8))') tetraData(n)%neighbors, sum(tetraData(n)%neighbors)
!         end do
    
    end subroutine assign_neighbors
    

    subroutine find_neighbors(vertIDS, list1, list2, tetraData)
    
        implicit none
        integer, intent(in)               :: vertIDs(:,:)
        integer, intent(inout)            :: list1(:), list2(:)
        type(tetraElement), intent(inout) :: tetraData(:)
        integer                           :: i, j, k, counter, face_i, face_j
        real                              :: t1, t2
        logical, dimension(4)             :: mask
   
        call cpu_time(t1)

        ! loop over elements in list1 and identify neighbors
        do i = 1,size(list1)
            
            ! count how many neighbors element i already has
            counter = count(tetraData(list1(i))%neighbors(:,1) > 0)
                                   
            ! search all elements in list2
            do j= 1,size(list2)
                
                ! look which vertices of element i are in the following j elements  
                mask = 0
                forall(k=1:4) mask(k) = any(vertIDs(k,list1(i)) == vertIDS(:,list2(j)))
                             
                ! 3 vertices are equal
                if (count(mask) == 3) then
                    
                    counter = counter + 1
                    
                    face_i = return_facenumber(mask)
                    
                    forall(k=1:4) mask(k) = any(vertIDs(k,list2(j)) == vertIDS(:,list1(i)))
                    face_j = return_facenumber(mask)
                    
                    tetraData(list1(i))%neighbors(face_i,1) = list2(j)
                    tetraData(list1(i))%neighbors(face_i,2) = face_j
                    tetraData(list2(j))%neighbors(face_j,1) = list1(i)
                    tetraData(list2(j))%neighbors(face_j,2) = face_i
                    
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
        
        call cpu_time(t2)
        write(*,'(a48,i8,a7,3es14.6)') "current search for neighbors with lists of size ", size(list1) + size(list2)," took: ",t2-t1
        write(*,*)
        
    end subroutine find_neighbors
    
    subroutine find_neighbors_single_list(vertIDS, list, tetraData)
    
        implicit none
        integer, intent(in)               :: vertIDs(:,:)
        integer, intent(inout)            :: list(:)
        type(tetraElement), intent(inout) :: tetraData(:)
        integer                           :: i, j, k, counter, face_i, face_j, n
        real                              :: t1, t2
        logical, dimension(4)             :: mask
   
        call cpu_time(t1)
        n = size(list)

        ! loop over elements in list1 and identify neighbors
        do i = 1,n-1
            
            ! if element has already 4 neighbors
            if (list(i) == 0) cycle
            
            ! count how many neighbors element i already has
            counter = count(tetraData(list(i))%neighbors(:,1) > 0)
                                   
            ! search all elements in list2
            do j= i+1,n
                
                ! look which vertices of element i are in the following j elements  
                mask = 0
                forall(k=1:4) mask(k) = any(vertIDs(k,list(i)) == vertIDS(:,list(j)))
                             
                ! 3 vertices are equal
                if (count(mask) == 3) then
                    
                    counter = counter + 1
                    
                    face_i = return_facenumber(mask)
                    
                    forall(k=1:4) mask(k) = any(vertIDs(k,list(j)) == vertIDS(:,list(i)))
                    face_j = return_facenumber(mask)
                    
                    tetraData(list(i))%neighbors(face_i,1) = list(j)
                    tetraData(list(i))%neighbors(face_i,2) = face_j
                    tetraData(list(j))%neighbors(face_j,1) = list(i)
                    tetraData(list(j))%neighbors(face_j,2) = face_i
                    
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
        
        call cpu_time(t2)
        write(*,'(a48,i8,a7,3es14.6)') "current search for neighbors within list of size ", n," took: ",t2-t1
        write(*,*)
        
    end subroutine find_neighbors_single_list
    
end module pre_process_data