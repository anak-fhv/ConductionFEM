! RayTracing modules
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 26.02.2014
! version: 0.2

module rt_constants

    implicit none
    integer, parameter   :: dp=selected_real_kind(p=14)
    real(dp), parameter  :: pi=3.14159265358979
end module rt_constants

module rt_types

    use rt_constants
    implicit none
    
    type :: tetraElement
        integer, dimension(4)    :: vertexIds=0   ! ID of vertices of the tetra element
        
        integer, dimension(4,2)  :: neighbors=0   ! information about of neighboring elements
        ! In above array the row index relates to the 4 faces of the element. 
        ! The 1st column contains the ID of the neighboring element or surface. In case of
        ! a surface the ID = ID_of_surface + Number_of_Tetraelements.
        ! The 2nd column contains the face of the neighboring element which is identical to
        ! the face given by the row index. In case a surface is the nieghbor the entry is -1.
        
        !real(dp), dimension(4,4) :: shape_funcs ! shape functions
        
        integer                  :: domain=0      ! to which domain the tetra belongs
    
    end type tetraElement
    
    type :: emissionSurface
        character(len=100)                   :: name ! name of emission surface
        real(dp), dimension(:), allocatable  :: area ! cumsum of area of the faces on the surface of emission
        integer, dimension(:,:), allocatable :: elemData ! (:,1) number of tetra-element      
                                                         ! (:,2) face of tetra-element which is on the surface
    end type
    
end module rt_types       

module math_funs
    
    use rt_constants
    implicit none
    
    contains
    
    ! put all zeros at the end of a list
    ! should contain only positive or zero entries
    subroutine putZerosAway(list)
        implicit none
        integer, intent(inout) :: list(:)
        integer                :: nmax, nzeros, value
        integer, dimension(1)  :: id
        
        nmax = size(list)
        nzeros = count(list == 0)
        if (nzeros == 0) return
 
        ! more of a debuuging feature, can be commented 
        if (minval(list) < 0) then
            write(*,*) "list should not contain negative entries!"
            write(*,*) list
            stop
        end if

        do
            id = minloc(list)  ! first index of smallest element 
            if (nmax-nzeros == id(1)-1) exit  ! all zero's are at the end
            
            ! move zero at the end, and shift remaining list 1 element to the front
            value = list(id(1))
            list(id(1):nmax-1) = list(id(1)+1:nmax) 
            list(nmax) = value
            
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
  
    ! wrapper for random numbers
    function myRandom(iflag)
        use ifport           ! used for standard random numbers           
        implicit none
        
        real(dp) :: myRandom
        integer  :: iflag
        
        if (iflag > 0) call srand(iflag)
        myRandom = drand(0)
        
    end function myRandom
    
end module math_funs                              


module rt_funcs

    use rt_types
    use rt_constants
    implicit none
    
    contains
    
    ! The function return_facenumber returns the ID of the face of a tetra
    ! given that a mask exist which indicates which vertices defining the face.        
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
    
    subroutine return_facevertIds(face, vertIds)
    
        implicit none
        integer, intent(in)                :: face
        integer, dimension(3), intent(out) :: vertIds
        
        select case(face)
            case(1)
                vertIds = [1,2,3]
            case(2)
                vertIds = [1,2,4]
            case(3)
                vertIds = [2,3,4]
            case(4)
                vertIds = [1,3,4]
            case default
                write(*,*) "Unknown face number!"
                stop
       end select
       
    end subroutine return_facevertIds    
    
   
end module rt_funcs


module pre_process_data

    use rt_funcs
    use math_funs
    
    implicit none
    
    contains
    
    subroutine read_mesh_data(file_name, emSurfNames, npart, tetraData, vertices, ems)
    
        implicit none
        character(len=*), intent(in)                                  :: file_name 
        character(len=*), intent(in)                                  :: emSurfNames(:)
        integer, intent(in)                                           :: npart
        type(tetraElement), dimension(:), allocatable, intent(out)    :: tetraData
        real(dp), dimension(:,:), allocatable, intent(out)            :: vertices 
        type(emissionSurface), dimension(:), allocatable, intent(out) :: ems
        logical                                                       :: l
        integer                                                       :: io_error, read_error, alloc_status, i, j, k, n
        integer                                                       :: nVertices, nTetra, nHexa, nPyr, nWedges, nDomain, nSurface
        integer, dimension(:,:), allocatable                          :: surfData
        integer, dimension(:), allocatable                            :: elemDomain
        real                                                          :: t1, t2
        integer                                                       :: nElem, remainElem
        character(len=100)                                            :: dName, surfName, line
        
    
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
        allocate(vertices(nVertices,3), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertex array")
        do n = 1,nVertices
            read(21,'(3(1x,3e14.6))', iostat=read_error) vertices(n,:) 
            call check_io_error(read_error,"reading vertex data")
        end do
        ! just some status messages
        write(*,'(a5,i7,a31,1x,i1,1x,i7)') "read ",nVertices," vertices into array of shape", shape(vertices) 
        write(*,*)
        
        ! allocate tetraData        
        allocate(tetraData(nTetra), stat=alloc_status)
        call check_alloc_error(alloc_status, "tetraData array")
        
        ! read tetra elements and assign vertices
        do n = 1,nTetra
            read(21,'(4(1x,i8))',iostat=read_error) tetraData(n)%vertexIds
            call check_io_error(read_error,"reading tetra elements")
        end do
        write(*,'(a5,i7,a)') "read ",nTetra," tetraeder elements" 
        write(*,*)
        
        ! read elements for each domain        
        do i = 1,nDomain
    
            ! read number of elements with domain and name of domain
            read(21,'(i8,1x,a8)', iostat=read_error) nElem, dName
            
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
        allocate(ems(size(emSurfNames)), stat=alloc_status)
        call check_alloc_error(alloc_status, "ems array")
        k = 0

        do i = 1,nSurface
    
            read(21,'(i8,1x,a)', iostat=read_error) nElem, surfName
        
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
            
            ! if the surface is a emission surface
            if (any(surfName == emSurfNames)) then
                k = k + 1
                allocate(ems(k)%area(nElem), stat=alloc_status)
                call check_alloc_error(alloc_status, "ems(k)%area array")
                allocate(ems(k)%elemData(nElem,2), stat=alloc_status)
                call check_alloc_error(alloc_status, "ems(k)%elemData array")
                
                ems(k)%name = surfName
                call CreateEmissionSurf(tetraData, vertices, surfData, ems(k))
                
            end if
            
            ! deallocate surfData
            deallocate(surfData, stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData, dealloaction") 
                
        end do
        
        ! assign neighbors
        call assign_neighbors(npart, tetraData)
        
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
    
    
    subroutine assign_neighbors(npart, tetraData)
    
        implicit none    
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
        
        ! loop over all partitions and perform searches
        do n = 1,npart
        
            ! search within same partition
            call find_neighbors_single_list(list(id_start(n):id_end(n)), tetraData)
!             write(*,*) "reduction after search in single list", count(list(id_start(n):id_end(n)) == 0)/real(id_end(n) - id_start(n) + 1)
!             write(*,*)
            
            
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
            do k = npart, n+1, -1
            
                call find_neighbors(list(id_start(n):id_end(n)), list(id_start(k):id_end(k)), tetraData)
!                 write(*,*) "reduction in high level list", count(list(id_start(n):id_end(n)) == 0)/real(id_end(n) - id_start(n) + 1)
!                 write(*,*) "reduction in low level list ", count(list(id_start(k):id_end(k)) == 0)/real(id_end(k) - id_start(k) + 1)
!                 write(*,*)
                
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
    
    end subroutine assign_neighbors
    

    subroutine find_neighbors(list1, list2, tetraData)
    
        implicit none
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
            
                ! if element has already 4 neighbors
                if (list2(j) == 0) cycle    
                
                ! look which vertices of element i are in the following j elements  
                mask = .false.
                forall(k=1:4) mask(k) = any(tetraData(list1(i))%vertexIds(k) == tetraData(list2(j))%vertexIds(:))
                             
                ! 3 vertices are equal
                if (count(mask) == 3) then
                    
                    counter = counter + 1
                    
                    face_i = return_facenumber(mask)
                    
                    ! determine the correct faces for element j
                    forall(k=1:4) mask(k) = any(tetraData(list2(j))%vertexIds(k) == tetraData(list1(i))%vertexIds(:))
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
        write(*,'(a49,i8,1x,i8,a7,3es14.6)') "current search for neighbors with lists of size ", size(list1),size(list2)," took: ",t2-t1
        write(*,*)
        
    end subroutine find_neighbors

        
    subroutine find_neighbors_single_list(list, tetraData)
    
        implicit none
        integer, intent(inout)            :: list(:)
        type(tetraElement), intent(inout) :: tetraData(:)
        integer                           :: i, j, k, counter, face_i, face_j, n
        real                              :: t1, t2
        logical, dimension(4)             :: mask
   
        call cpu_time(t1)
        n = size(list)        

        
        ! loop over elements in list and identify neighbors
        do i = 1,n-1
                 
                                    
            ! if element has already 4 neighbors
            if (list(i) == 0) cycle
            
            ! count how many neighbors element i already has
!             write(*,*) i, n-1
!             write(*,*)
            counter = count(tetraData(list(i))%neighbors(:,1) > 0)
                                   
            ! search all elements below
            do j= i+1,n
            
                ! if element has already 4 neighbors
                if (list(j) == 0) cycle
                
                ! look which vertices of element i are in the following j elements  
                mask = .false.
                forall(k=1:4) mask(k) = any(tetraData(list(i))%vertexIds(k) == tetraData(list(j))%vertexIds(:))
                             
                ! 3 vertices are equal
                if (count(mask) == 3) then
                    
                    counter = counter + 1
                    
                    face_i = return_facenumber(mask)
                    
                    forall(k=1:4) mask(k) = any(tetraData(list(j))%vertexIds(k) == tetraData(list(i))%vertexIds(:))
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
        write(*,'(a47,i12,a7,3es14.6)') "current search for neighbors with list of size ", n," took: ",t2-t1
        write(*,*)
        
    end subroutine find_neighbors_single_list

        
    subroutine CreateEmissionSurf(tetras, vertices, surfData, ems)

        implicit none
        type(tetraElement), intent(in)            :: tetras(:)
        real(dp), intent(in)                      :: vertices(:,:)   
        integer, intent(in)                       :: surfData(:,:)
        type(emissionSurface), intent(inout)      :: ems        
        integer                                   :: n,i
        integer, dimension(3)                     :: vertIds
        real(dp), dimension(3)                    :: p1, p2, p3 
        
        do i =1,size(surfData,1)
            
                ! face and the respective vertices which lie on the surface of emission
                call return_facevertIds(surfData(i,2),vertIds)           
                p1 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(1)),:)
                p2 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(2)),:)
                p3 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(3)),:)
            
                ! get area normal vectors
                ems%ElemData(i,1) = surfData(i,1)
                ems%ElemData(i,2) = surfData(i,2)
                if (i == 1) then
                    ems%area(1) = 0.5_dp*norm(cross(p2-p1,p3-p1))
                else 
                    ems%area(i) = 0.5_dp*norm(cross(p2-p1,p3-p1)) + ems%area(i-1)
                end if 
        end do
        
        ems%area = ems%area/ems%area(size(surfData,1))
        
    end subroutine CreateEmissionSurf
    

end module pre_process_data