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
        integer, dimension(4)    :: vertexIds=0   ! ID ofvertex vertices of the tetra element
        
        integer, dimension(4,2)  :: neighbors=0   ! information about of neighboring elements
        ! In above array the row index relates to the 4 faces of the element. 
        ! The 1st column contains the ID of the neighboring element or surface. In case of
        ! a surface the ID = ID_of_surface + Number_of_Tetraelements.
        ! The 2nd column contains the face of the neighboring element which is identical to
        ! the face given by the row index. In case a surface is the nieghbor the entry is -1.
        
        !real(dp), dimension(4,4) :: shape_funcs ! shape functions
        
        integer                  :: domain=0      ! to which domain the tetra belongs
        integer                  :: nAbsorbed=0   ! number of absorbed elements
    
    end type tetraElement
    
    type :: emissionSurface
        character(len=100)                   :: name ! name of emission surface
        real(dp), dimension(:), allocatable  :: area ! cumsum of area of the faces on the surface of emission
        integer, dimension(:,:), allocatable :: elemData ! (:,1) number of tetra-element      
                                                         ! (:,2) face of tetra-element which is on the surface
    end type
    
    type :: rayContainer  ! contains information on the traced ray
        real(dp), dimension(3) :: point, direction ! current point and its direction
        integer                :: tetraID          ! current tetraeder indes
        integer                :: faceID           ! current face index
        real(dp)               :: length=0.0_dp    ! distance travelled
    end type
    
end module rt_types       


module math_funs
    
    use rt_constants
    implicit none
    
    contains
    
    ! put all zeros at the end of a list
    ! list should contain only positive or zero entries
    subroutine putZerosAway(list)

        integer, intent(inout) :: list(:)
        integer                :: nmax, nzeros, value
        integer, dimension(1)  :: id
        
        nmax = size(list)
        nzeros = count(list == 0)
        if (nzeros == 0) return
 
        ! more of a debuuging feature, can be commented size
        if (minval(list) < 0) then
            write(*,*) "list should not contain negative entries!"
            write(*,*) list
            stop
        end if

        do
            id = minloc(list)                 ! first index of smallest element 
            if (nmax-nzeros == id(1)-1) exit  ! all zero's are at the end exit loop
            
            ! move zero at the end, and shift remaining list about 1 element to the front
            value = list(id(1))
            list(id(1):nmax-1) = list(id(1)+1:nmax) 
            list(nmax) = value
            
        end do
        
    end subroutine putZerosAway

    ! determine a rotation matrix given an axis of rotation (rotAxis),
    ! an angle of rotation (rotAngle), and a vector perpendicular (vecPerp) to it
    subroutine RotationMatrix(rotAxis, vecPerp, rotAngle, rotMatrix)
    
         real(dp), dimension(3), intent(in)    :: rotAxis, vecPerp
         real(dp)                              :: rotAngle
         real(dp), dimension(3,3), intent(out) :: rotMatrix
         real(dp), dimension(3,3)              :: A,M
         
         ! cartesian basis
         A = reshape([rotAxis, vecPerp, cross(rotAxis, vecPerp)],[3,3])
         ! general rotation matrix
         M = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(rotAngle), sin(rotAngle), 0.0_dp, -sin(rotAngle), cos(rotAngle)],[3,3])
         !rotation matrix for cartesian basis
         rotMatrix = matmul(matmul(A,M),transpose(A))
         
    end subroutine RotationMatrix
    
    ! cross product in 3d
    function cross(v1,v2)

        real(dp), dimension(3) :: cross
        real(dp), dimension(3), intent(in)  :: v1,v2
 
        cross(1) = v1(2)*v2(3) - v1(3)*v2(2) 
        cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
 
    end function cross
  
    ! norm 
    function norm(v)

        real(dp)             :: norm
        real(dp), intent(in) :: v(:)
  
        norm = sqrt(dot_product(v,v))

   end function norm
  
    ! wrapper for random numbers
    function myRandom(iflag)
        use ifport           ! use intel random numbers           
        
        real(dp) :: myRandom
        integer  :: iflag
        
        if (iflag > 0) call srand(iflag) ! seeds random number
        myRandom = drand(0)
        
    end function myRandom
    
end module math_funs                              


module rt_funcs

    use rt_types
    use math_funs
    implicit none
    
    contains
    
    ! The function return_facenumber returns the ID of the face of a tetra
    ! given that a mask exist which indicates which vertices define the face.        
    integer function return_facenumber(mask)
            
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
    
    ! returns the indices of tetraeder points (between 1 and 4) for a given face 
    subroutine return_facevertIds(face, vertIds)
    
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
    
    ! test whether a point is inside a triangle or not
    function PointInside(p1,p2,p3,tp)
        
        real(dp), dimension(3), intent(in) :: p1,p2,p3,tp
        logical                            :: PointInside
        real(dp), dimension(3)             :: na,nx,ny,nz,bc
        real(dp)                           :: dpna
        
        na = cross(p2-p1,p3-p1)
        nx = cross(p3-p2,tp-p2)
        ny = cross(p1-p3,tp-p3)
        nz = cross(p2-p1,tp-p1)
        dpna = dot_product(na,na)
        bc(1) = dot_product(na,nx)/dpna
        bc(2) = dot_product(na,ny)/dpna
        bc(3) = dot_product(na,nz)/dpna
        PointInside = (all(bc > 0e0_dp) .and. (abs(sum(bc) - 1) <= 1e-14_dp))
        
!         write(*,*) bc(1)
!         write(*,*) bc(2)
!         write(*,*) bc(3)
!         write(*,*) 1-sum(bc)
    
    end function
    
   
end module rt_funcs


module helper_functions

    implicit none
    
    contains
    
        ! checks for input-output errors
    subroutine check_io_error(stat, message, unitnumber)
    
        integer, intent(in)                :: stat, unitnumber
        character(*), intent(in), optional :: message
        
        if ((stat > 0) .and. present(message)) then
            write(*,'(A)') "io-error: "//message
            close(unit=unitnumber)
            stop
        else if (stat > 0) then
            write(*,'(A)') "error in input-output/reading/writing procedure!"
            close(unit=unitnumber)
            stop
        end if
        
        return
    end subroutine check_io_error

    ! checks for alloaction errors    
    subroutine check_alloc_error(stat, message)
    
        integer, intent(in)                :: stat
        character(*), intent(in), optional :: message
        
        if ((stat /= 0) .and. present(message)) then
            write(*,'(A)') "allocation-error: "//message
            stop
        else if (stat > 0) then
            write(*,'(A)') "error in allocation procedure!"
            stop
        end if
        
        return
    end subroutine check_alloc_error
    
end module helper_functions


module pre_process_data

    use rt_funcs
    use math_funs
    use helper_functions
    
    implicit none
    
    contains
    
    ! main routine which perform reading of input data and calls routines for
    ! creating emission surfaces and populates tetraElement-type
    subroutine read_mesh_data(file_name, emSurfNames, npart, tetraData, vertices, ems)

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
        write(*,*)
        write(*,*) "Start Reading Input File"
        write(*,*) "==========================="
        write(*,*)
    
        ! open file
        open(unit=81, file=file_name, status='old', action='read', iostat=io_error)       
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
            read(81,'(3(1x,3e14.6))', iostat=read_error) vertices(n,:) 
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
        allocate(ems(size(emSurfNames)), stat=alloc_status)
        call check_alloc_error(alloc_status, "ems array")
        k = 0

        do i = 1,nSurface
    
            ! read number of elements with domain and name of domain
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
            ! surface index starts at number of (tetra elements + 1)
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
                
                ! create emission surface
                ems(k)%name = surfName
                call CreateEmissionSurf(tetraData, vertices, surfData, ems(k))
                
            end if
            
            ! deallocate surfData
            deallocate(surfData, stat=alloc_status)
            call check_alloc_error(alloc_status, "surfData, dealloaction") 
                
        end do
        
        ! close file
        close(unit=81)
        
        ! assign neighbors
        write(*,*) "Creating Element Connections"
        write(*,*) "============================="
        write(*,*) 
        call assign_neighbors(npart, tetraData)
        
        ! check for double entries (could be commented)
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
        
        ! report timing
        call cpu_time(t2)
        write(*,*)
        write(*,*) "Finished Pre-Processing"
        write(*,*) "=============================" 
        write(*,*) "pre-processing took (in Sec.): ", t2-t1
    
    end subroutine read_mesh_data


    ! main routine for assigning neighbors
    subroutine assign_neighbors(npart, tetraData)
       
        integer, intent(in)                            :: npart
        type(tetraElement), intent(inout)              :: tetraData(:)
        integer                                        :: alloc_status, n, k, r
        integer, dimension(:), allocatable             :: list, id_start, id_end
        integer, dimension(1)                          :: idFirstZero, maxElem
        
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
        end do
        
        write(*,*) "number of partitions ", npart
        write(*,*) "number of valid elements ", maxElem
!         write(*,*) id_start
!         write(*,*) id_end
        
        ! loop over all partitions and perform searches
        do n = 1,npart
        
            write(*,*) "parsing partition number", n 
            ! search within same partition
            call find_neighbors_single_list(list(id_start(n):id_end(n)), tetraData)            
            
            ! filter list of partition n by putting zeros at the end of the list
            call filter_list(list(id_start(n):id_end(n)),id_start(n), id_end(n))
            if (id_end(n) == 0) exit
            
            ! search in all partitions with larger index number
            do k = npart, n+1, -1
            
                call find_neighbors(list(id_start(n):id_end(n)), list(id_start(k):id_end(k)), tetraData)
                
                ! filter list of partition n
                call filter_list(list(id_start(n):id_end(n)),id_start(n), id_end(n))
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
    subroutine find_neighbors(list1, list2, tetraData)

        integer, intent(inout)            :: list1(:), list2(:)
        type(tetraElement), intent(inout) :: tetraData(:)
        integer                           :: i, j, k, counter, face_i, face_j
        logical, dimension(4)             :: mask
!         real                              :: t1, t2
   
!         call cpu_time(t1)

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
                    call add_neighbor(list1(i), list2(j), face_i, tetraData)                     
                    
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
        
!         call cpu_time(t2)
!         write(*,'(a49,i8,1x,i8,a7,3es14.6)') "current search for neighbors with lists of size ", size(list1),size(list2)," took: ",t2-t1
!         write(*,*)
        
    end subroutine find_neighbors

    ! search for neighbors in single lists
    ! search routine:
    ! for each element in the list a neighbor is looked for in all following elements    
    subroutine find_neighbors_single_list(list, tetraData)

        integer, intent(inout)            :: list(:)
        type(tetraElement), intent(inout) :: tetraData(:)
        integer                           :: i, j, k, counter, face_i, face_j, n
        logical, dimension(4)             :: mask
!         real                              :: t1, t2
   
!         call cpu_time(t1)
        n = size(list)        

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
                    call add_neighbor(list(i), list(j), face_i, tetraData)
                    
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
        
!         call cpu_time(t2)
!         write(*,'(a49,i12,a7,3es14.6)') "current search for neighbors with list of size ", n," took: ",t2-t1
!         write(*,*)
        
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
    subroutine add_neighbor(tetra1, tetra2, face1, tetraData)
    
        integer, intent(in)               :: tetra1, tetra2, face1
        type(tetraElement), intent(inout) :: tetraData(:)
        logical, dimension(4)             :: mask
        integer                           :: k, face2
        
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
    subroutine CreateEmissionSurf(tetras, vertices, surfData, ems)

        type(tetraElement), intent(in)            :: tetras(:)
        real(dp), intent(in)                      :: vertices(:,:)   
        integer, intent(in)                       :: surfData(:,:)
        type(emissionSurface), intent(inout)      :: ems        
        integer                                   :: n,i
        integer, dimension(3)                     :: vertIds
        real(dp), dimension(3)                    :: p1, p2, p3 
        
        do i =1,size(surfData,1)
            
                ! index of tetraeder and the index for the face on the surface 
                ems%elemData(i,1) = surfData(i,1)
                ems%elemData(i,2) = surfData(i,2)
                
                ! determine vertices for calculating the face area
                call return_facevertIds(surfData(i,2),vertIds)           
                p1 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(1)),:)
                p2 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(2)),:)
                p3 = vertices(tetras(surfData(i,1))%vertexIds(vertIds(3)),:)                
                if (i == 1) then
                    ems%area(1) = 0.5_dp*norm(cross(p2-p1,p3-p1))
                else 
                    ems%area(i) = 0.5_dp*norm(cross(p2-p1,p3-p1)) + ems%area(i-1)
                end if 
                
        end do
        
        ! normalize area
        ems%area = ems%area/ems%area(size(surfData,1))
        
    end subroutine CreateEmissionSurf
    
    ! write out data from pre-processing
    subroutine WriteData(tetraData, vertices, emSurf, file_name)
    
        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(emissionSurface), intent(in)   :: emSurf(:)
        character(len=*), intent(in)        :: file_name
        character(len=100)                  :: fname
        integer                             :: io_error, n, k, write_error
        
        ! create file from exisiting file-name
        fname = "../obj/"//file_name(1:index(file_name,".msh")-1)//".dat"
        open(unit=82, file=fname, status='replace', action='write', iostat=io_error)       
        call check_io_error(io_error,"opening file for write-out",82)
        
        ! write values for number of various elements
        write(82,'(3(1x,i9))') size(vertices,1), size(tetraData), size(emSurf)
        
        ! write vertices
        do n = 1,size(vertices,1)
            write(82,'(3(1x,3e14.6))', iostat=write_error) vertices(n,:) 
            call check_io_error(write_error,"writing vertex data",82)
        end do
        
        
        ! write tetraData (vertices and domain)
        do n = 1,size(tetraData)
            write(82,'(5(1x,i8))',iostat=write_error) tetraData(n)%vertexIds, tetraData(n)%domain 
            call check_io_error(write_error,"writing tetra data",82)
        end do
        
        ! write tetraData (connections)
        do n = 1, size(tetraData)
            write(82,'(1x,4(i8,1x,i2))',iostat=write_error) transpose(tetraData(n)%neighbors) 
            call check_io_error(write_error,"writing tetra data part 2",82)
        end do
        
        ! write emission surface information
        do n = 1, size(emSurf)
            write(82,'(1x,i8,1x,a)',iostat=write_error) size(emSurf(n)%area), trim(emSurf(n)%name)
            call check_io_error(write_error,"writing emission surface info",82)
            
            do k=1,size(emSurf(n)%area)
                write(82,'(1x,i8,1x,i2,2x,3e14.6)',iostat=write_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2), emSurf(n)%area(k)
                call check_io_error(write_error,"writing emission surface data",82)
            end do
        end do    
        
        ! close file
        close(unit=82)
  
    end subroutine WriteData
   
    ! read data from exisiting pre-processing file
    subroutine ReadData(tetraData, vertices, emSurf,file_name)
    
        type(tetraElement), allocatable, dimension(:), intent(out)    :: tetraData
        real(dp), allocatable, dimension(:,:), intent(out)            :: vertices
        type(emissionSurface), allocatable, dimension(:), intent(out) :: emSurf
        character(len=*), intent(in)                                  :: file_name
        integer :: io_error, n, k, read_error, nVertices, nTetra, nSurf, alloc_status, nElem                                        
        logical :: l
        
        ! check if file exists
        inquire(file="../obj/"//file_name, exist=l)
        if (l .eqv. .false.) then
            write (*,*) file_name//" does not exist!"
            stop
        end if
        
        ! open file
        open(unit=83, file="../obj/"//file_name, status='old', action='read', iostat=io_error)       
        call check_io_error(io_error,"opening file for reading pre-processed data",83)
        
        ! read information about element sizes
        read(83,'(3(1x,i9))') nVertices, nTetra, nSurf
        
        ! allocate vertices and read vertices
        allocate(vertices(nVertices,3), stat=alloc_status)
        call check_alloc_error(alloc_status, "vertex array")
        do n = 1,nVertices
            read(83,'(3(1x,3e14.6))', iostat=read_error) vertices(n,:) 
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
            read(83,'(1x,4(i8,1x,i2))',iostat=read_error) (tetraData(n)%neighbors(k,1), tetraData(n)%neighbors(k,2), k=1,4)
            call check_io_error(read_error,"reading tetra data part 2",83)
        end do
        
        ! allocate emission surface and read information
        allocate(emSurf(nSurf), stat=alloc_status)
        call check_alloc_error(alloc_status, "emSurf array")
        do n = 1, nSurf
            read(83,'(1x,i8,1x,a)',iostat=read_error) nElem, emSurf(n)%name
            call check_io_error(read_error,"writing emission surface info",83)
            
            allocate(emSurf(n)%area(nElem), stat=alloc_status)
            call check_alloc_error(alloc_status, "emSurf%area array")
            allocate(emSurf(n)%elemData(nElem,2), stat=alloc_status)
            call check_alloc_error(alloc_status, "emSurf%elemData array")
            
            do k=1,nElem
                read(83,'(1x,i8,1x,i2,2x,3e14.6)',iostat=read_error) emSurf(n)%elemData(k,1), emSurf(n)%elemData(k,2), emSurf(n)%area(k)
                call check_io_error(read_error,"writing emission surface data",83)
            end do
        end do    
        
        ! close file
        close(unit=83)
  
   end subroutine ReadData

end module pre_process_data


module tracing

    use rt_funcs
    use math_funs
    use helper_functions
    
    implicit none
    
    contains
    
    ! given an emission surface, the following subroutine performs the following steps:
    ! 1. select a face on the emission surface by roulette wheel selection
    ! 2. select a point within the triangular face based on 2 random numbers 
    !    (uses triangle and uniform distribution)
    ! 3. select a random direction    
    subroutine CreateRay(tetraData, vertices, ems, ray)

        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(emissionSurface), intent(in)   :: ems
        type(rayContainer), intent(out)     :: ray
        integer                             :: i
        integer, dimension(3)               :: vertIDs 
        real(dp)                            :: psi, theta, b, c, d
        real(dp), dimension(3)              :: p1, p2, p3, dir1, dir2, dir21, ds1, ndir
        real(dp), dimension(3,3)            :: M  

        ! get random tetraeder on the emission surface
        psi = myRandom(0) ! initialize random generator
        
        ! decision whether to start at the beginning or end of the list (just for speed)
        if (psi > 0.5) then
            do i = size(ems%area)-1,1,-1
                if (ems%area(i) < psi) exit
            end do
            ray%tetraID = ems%elemData(i+1,1)
            ray%faceID = ems%elemData(i+1,2)
        else           
            do i = 1,size(ems%area)
                if (ems%area(i) > psi) exit
            end do
            ray%tetraID = ems%elemData(i,1)
            ray%faceID = ems%elemData(i,2)
        end if
    
        ! get vertices for the selected face on the surface
        call return_facevertIds(ray%faceID,vertIDs)  
        p1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
        p2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:)
        p3 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:)
     
        ! find longest side of triangle
        if ((norm(p2-p1) >= norm(p3-p1)) .and. (norm(p2-p1) >= norm(p3-p2))) then
            dir1 = p2-p1
            dir2 = p3-p1
            ray%point = p1     
!             write(*,*) "p1"
!             write(*,*)   
        else if ((norm(p1-p3) >= norm(p2-p1)) .and. (norm(p1-p3) >= norm(p3-p2))) then
            dir1 = p1-p3
            dir2 = p2-p3 
            ray%point = p3
!             write(*,*) "p3"
!             write(*,*)  
        else
            dir1 = p3-p2
            dir2 = p1-p2
            ray%point = p2
!             write(*,*) "p2"
!             write(*,*)
        end if
        
        ! projection of dir2 onto dir1
        dir21 = dot_product(dir1,dir2)/dot_product(dir1,dir1) * dir1
    
        ! choose random numbers based on triangle distribution
        psi = myRandom(0)
        theta = myRandom(0)

        b = norm(dir1)     ! length of dir1
        c = norm(dir21)    ! length of projection dir21
        ds1 = dir2 - dir21 ! perpendicular to d1 and going through vertex of triangle which is not on d1

        if (psi <= c/b) then
            ray%point = ray%point + sqrt(psi*b*c)*dir1/b                       ! random number along longest side (triangle distribution)
            d = sqrt(psi*b*c)*tan(acos(dot_product(dir1,dir2)/(b*norm(dir2)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d*ds1/norm(ds1)                      ! random number along perpendicular direction (uniform distribution)     
        else
            ray%point = ray%point + (b - sqrt(b*(b-c)*(1-psi)))*dir1/b                               ! random number along longest side (triangle distribution)
            d = sqrt(b*(b-c)*(1-psi)) *tan(acos(dot_product(-dir1,(dir2-dir1))/(b*norm(dir2-dir1)))) ! length of perpendicular vector at new rayOrigin
            ray%point = ray%point + theta*d*ds1/norm(ds1)                                            ! random number along perpendicular direction (uniform distribution)        
        end if
        
        ! test whether point is indeed in triangle
        if (PointInside(p1,p2,p3,ray%point) .eqv. .false.) then
            write(*,*) "Emission Point not in chosen face!"
            stop
        end if
    
        ! choose ray direction
        ndir = cross(dir1,dir2) ! normal vector of triangle pointing outwards
        ndir = ndir/norm(ndir)
        
        ! check whether normal vector points inwards or outwards
        do i = 1,4
            if (any(i == vertIDs) .eqv. .false.) exit
        end do        
        if (dot_product(ndir, vertices(tetraData(ray%tetraID)%vertexIds(i),:) - ray%point) < 0) ndir = -ndir

        ! first rotate about about a vector in the plane
        dir1 = dir1/b
        call RotationMatrix(dir1, ndir, asin(sqrt(myRandom(0))), M)
        ray%direction = matmul(M,ndir)
        
        ! second rotation with normal vector as rotation axis
        call RotationMatrix(ndir, dir1, asin(sqrt(myRandom(0))), M)
        ray%direction = matmul(M,ray%direction)
        
!         write(*,*)
!         write(*,'(a12,3(1x,3e14.6))') "origin: ", ray%point
!         write(*,'(a12,3(1x,3e14.6))') "direction: ", ray%direction
!         write(*,*)         
        
    end subroutine CreateRay
    
    ! main raytracing routine
    subroutine TraceRay(tetraData, vertices, ray)
    
        type(tetraElement), intent(inout)   :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        real(dp)                            :: kappa, sigma ! absorption and scattering coefficients (must be provided from somewhere)
        real(dp)                            :: lAbs, lScat
        integer, dimension(3)               :: vertIDs 
        real(dp), dimension(3)              :: rp, v1, v2
        logical                             :: test
        integer                             :: counter=0, io_error, write_error
        
        ! calculate length of initial ray
        kappa = 1.0_dp ! simple assumption so far
        sigma = 1.0_dp ! simple assumption so far
        
        lAbs = 1.0_dp/kappa*log(1/myRandom(0))
        lScat = 1.0_dp/sigma*log(1/myRandom(0))
        
        open(unit=84, file="../obj/sphere-absorbed.out", action='write', iostat=io_error, position='append')       
        call check_io_error(io_error,"opening file for results 1",84)
                
        do
            if (ray%length >= lAbs) then
                write(*,*) "ray is absorbed"
                tetraData(ray%tetraID)%nAbsorbed = tetraData(ray%tetraID)%nAbsorbed+1
                write(84,'(1x,i8,1x,3(3e14.6,1x),3e14.6)',iostat=write_error) ray%tetraID, ray%point, myRandom(0)
                call check_io_error(write_error,"writing results 1",84)
                counter = counter + 1
                exit
            end if
            
            call FindNextTetra(vertices,tetraData,ray)
            
!             write(*,*) ray%faceID
!             write(*,*) ray%tetraID
!             write(*,*) ray%length
!             write(*,*) ray%point
!             write(*,*)
        
            if (ray%faceID == -1) then
                write(*,*) "point on enclosure"
                
                exit
            end if
        
            ! just some test
            call return_facevertIds(ray%faceID, vertIDs)
            rp = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
            v1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:)   
            v2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:)
            test = PointInside(rp,v1,v2,ray%point)
            
            if (test .eqv. .false.) then
                write(*,*) "POINT NOT ON TETRA FACE!!!"
                stop
            end if
            
        end do
        
        close(unit=84, iostat=io_error)
        call check_io_error(io_error,"opening file for results 1",84)
        
    end subroutine TraceRay

    ! find next tetra Element along ray direction
    subroutine FindNextTetra(vertices,tetraData,ray)
        
        type(tetraElement), intent(in)      :: tetraData(:)
        real(dp), intent(in)                :: vertices(:,:)
        type(rayContainer), intent(inout)   :: ray
        real(dp) :: alpha, temp
        integer  :: f, newFace
        real(dp), dimension(3) :: rp, v1, v2, nsf
        integer, dimension(3)  :: vertIDs
        
        alpha = 100.0_dp  ! some large initial value for alpha
        do f = 1,4
        
            if (f == ray%faceID) cycle  ! is face where ray is emitted
            
            ! 1 vertex and 2 vectors of current face
            call return_facevertIds(f,vertIDs)
            rp = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(1)),:)
            v1 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(2)),:) - rp  
            v2 = vertices(tetraData(ray%tetraID)%vertexIds(vertIDs(3)),:) - rp
            
            ! normal vector of current face  
            nsf = cross(v1,v2)
            nsf = nsf/norm(nsf) 

            ! calculate length until intersection
            temp = (dot_product(nsf,rp) - dot_product(nsf,ray%point))/dot_product(nsf, ray%direction)
            if (temp < 0) cycle ! length can not be negative
            
            ! check if temp is the shortest length so far
            ! if true, current face is the face the ray will intersect
            if (temp < alpha) then
                alpha = temp
                newFace = f
            end if
            
        end do
        
        ! update ray container
        ray%length = ray%length+norm(alpha*ray%direction)
        ray%point = ray%point + alpha*ray%direction
        ray%faceID = tetraData(ray%tetraID)%neighbors(newFace,2)  
        ray%tetraID = tetraData(ray%tetraID)%neighbors(newFace,1)     
    
    end subroutine
        
    ! calculate refraction with Snell's law
    subroutine SnellsLaw(incident, n1, n2, nsf, reflect, refract)
        
        real(dp), intent(in)                  :: n1, n2   ! refractive indices
        real(dp), dimension(3), intent(in)    :: incident ! incident direction 
        real(dp), dimension(3), intent(inout) :: nsf      ! normal vector
        real(dp), dimension(3), intent(out)   :: refract, reflect ! refraction and reflection direction
        real(dp) :: cosAngle, ratio
        
        cosAngle = -dot_product(incident,nsf)
        if (cosAngle < 0) then
            cosAngle = -cosAngle
            nsf = -nsf
        end if
        
        reflect = incident + 2*cosAngle*nsf
        
        ratio = n1/n2
        refract = ratio*incident + ratio*cosAngle*nsf - sqrt(1 - ratio**2*(1-cosAngle**2))*nsf
   end subroutine SnellsLaw
        
end module tracing        
    
    
    