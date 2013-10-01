module numericalmethods
    implicit none
    contains

    subroutine solvedeterminant(inputmatrix,determinant)
        double precision, dimension(:,:) :: inputmatrix
        double precision, dimension(:,:), allocatable :: funcmat
        double precision :: determinant
        integer :: i,inputsize
        inputsize = size(inputmatrix,1)
        allocate(funcmat(inputsize,inputsize))
        funcmat = inputmatrix
        determinant = gaussiandeterminant(funcmat)
        deallocate(funcmat)
    end subroutine solvedeterminant

    function simpleinverse(A) result(B)
        double precision,dimension(:,:) :: A
        double precision,dimension(:,:),allocatable :: B,Cof
        double precision :: dm
        integer :: i,j,k,n

        n = size(A,1)
        allocate(B(n,n))

        do i=1,n
            do j=1,n
                Cof(i,j) = (-1**(i+j))*cofactor(A,i,j)
            end do
        end do
    end function simpleinverse

    function cofactor(A,i,j) result(co)
        double precision,dimension(:,:) :: A
        double precision,dimension(:,:), allocatable :: B
        double precision :: co
        integer :: i,j,k,n

        n = size(A)
        co = 0
        allocate(B(n,n))


    end function cofactor

    function getminor(inputmatrix,i,j) result(minor)
        integer :: i,j
        integer :: m,p,q
        double precision, dimension(:,:) :: inputmatrix
        double precision, dimension(:,:), allocatable :: temp1,minor
        m = size(inputmatrix,1)
        allocate(temp1(m-1,m))
        allocate(minor(m-1,m-1))
        q = 0
        do p=1,m
            if(p == i) then
                continue
            else
                q = q+1
                temp1(q,:) = inputmatrix(p,:)
            end if
        end do
        q = 0
        do p=1,m
            if(p == j) then
                continue
            else
                q = q+1
                minor(:,q) = temp1(:,p)
            end if
        end do

        return
    end function getminor

    function gaussiandeterminant(A) result(dm)
        double precision, dimension(:,:) :: A
        double precision :: dm
        integer :: i,j,n,p
        integer,dimension(1) :: maxl
        n = size(A,1)
        do i=1,n-1
            maxl = maxloc(abs(A(i:n,i)))
            p = maxl(1) + i -1
            if(p /= i) then
                A = pivot(A,i,p)
            end if
            do j=i+1,n
                A(j,:) = A(j,:) - (A(j,i)/A(i,i))*A(i,:)
            end do
        end do
        dm = 1.0
        do i=1,n
            dm = dm*A(i,i)
        end do
    end function gaussiandeterminant

    function sortedbinarysearch(A,b) result(indexof)
        integer,dimension(:) :: A
        integer :: b,indexof,first,last,mid
        first = 1
        last = size(A,1)
        indexof = 0
        do while(last-first >= 1)
            mid = (last+first)/2
            print *, 'Last: ', last, 'First: ', first, 'Mid: ', mid
            if(b == A(mid)) then
                indexof = mid
                return
            elseif(b < A(mid)) then
                last = mid -1
            else
                first = mid + 1
            end if
        end do
        return
    end function sortedbinarysearch

    subroutine gaussjordaninverse(inputmatrix,inversematrix)
        double precision,dimension(:,:) :: inputmatrix, inversematrix
        double precision,dimension(:,:),allocatable :: funcmat
        integer :: inputsize,i,j,k

        inputsize = size(inputmatrix,1)
        allocate(funcmat(inputsize,inputsize))
        funcmat = inputmatrix
        inversematrix = partialpivotedinverse(funcmat,inputsize)
    end subroutine gaussjordaninverse

    function nonpivotinverse(A,n) result(B)
        double precision,dimension(:,:) :: A
        double precision,dimension(:,:),allocatable :: B, C
        double precision :: factor
        integer :: i,j,k,n,p
        integer,dimension(1) :: maxl

        allocate(B(n,n))
        allocate(C(n,2*n))
        B = eye(n)
        C(1:n,1:n) = A
        C(1:n,n+1:2*n) = B
        do i=1,n-1
            do j=i+1,n
                factor = C(j,i)/C(i,i)
                C(j,:) = C(j,:) - factor*C(i,:)
                C(i,:) = C(i,:)/C(i,i)
            end do
        end do
        C(n,:) = C(n,:)/C(n,n)
        do i=n,2,-1
            do j=i-1,1,-1
                factor = C(j,i)/C(i,i)
                C(j,:) = C(j,:) - factor*C(i,:)
            end do
        end do
        B = C(1:n,n+1:2*n)
        deallocate(C)
    end function nonpivotinverse

    function partialpivotedinverse(A,n) result(B)
        double precision,dimension(:,:) :: A
        double precision,dimension(:,:),allocatable :: B,C
        double precision :: factor
        integer :: i,j,k,n,p
        integer,dimension(1) :: maxl

        allocate(B(n,n))
        allocate(C(n,2*n))
        B = eye(n)
        C(1:n,1:n) = A
        C(1:n,n+1:2*n) = B
        do i=1,n-1
            maxl = maxloc(abs(C(i:n,i)))
            p = maxl(1) + i -1
            if(p /= i) then
                C = pivot(C,i,p)
            end if
            do j=i+1,n
                factor = C(j,i)/C(i,i)
                C(j,:) = C(j,:) - factor*C(i,:)
                C(i,:) = C(i,:)/C(i,i)
            end do
        end do
        C(n,:) = C(n,:)/C(n,n)
        do i=n,2,-1
            do j=i-1,1,-1
                factor = C(j,i)/C(i,i)
                C(j,:) = C(j,:) - factor*C(i,:)
            end do
        end do
        B = C(:,n+1:2*n)
        deallocate(C)
    end function partialpivotedinverse

    function pivot(A,i,j) result(M)
        double precision,dimension(:,:) :: A
        double precision,dimension(size(A,1),size(A,2)):: M
        double precision,dimension(:),allocatable :: temprow
        integer :: i,j,rows,cols

        rows = size(A,1)
        cols = size(A,2)
        allocate(temprow(cols))
        M = A
        temprow = M(j,:)
        M(j,:) = M(i,:)
        M(i,:) = temprow
        deallocate(temprow)
        return
    end function pivot

    function eye(n) result(U)
        integer :: n,i
        double precision,dimension(:,:),allocatable :: U
        allocate(U(n,n))
        U = 0.0d0
        do i=1,n
            U(i,i) = 1
        end do
        return
    end function eye

    function ludecomp(A) result(B)
        double precision,dimension(:,:) :: A
        double precision,dimension(size(A,1),size(A,1)) :: B
        integer :: i,j,k,n
        double precision :: factor,innersum

        n = size(A,1)
        B = eye(n)
        do i=1,n-1
            do j=i+1,n
                factor = A(j,i)/A(i,i)
                A(j,i:n) = A(j,i:n) - factor*A(i,i:n)
                A(j,i) = factor
            end do
        end do

        do i=1,n
            do j=1,n
                if(j==1) then
                    B(j,i) = A(j,j)*B(j,i)/A(j,j)
                else
                    innersum = 0.0d0
                    do k=1,j-1
                        innersum = innersum + A(j,k)*B(k,i)
                    end do
                    B(j,i) = (A(j,j)/A(j,j))*(B(j,i) - innersum)
                end if
            end do
            do j=n,1,-1
                if(j==n) then
                    B(j,i) = B(j,i)/A(j,j)
                else
                    innersum = 0.0d0
                    do k=n,j+1,-1
                        innersum = innersum + A(j,k)*B(k,i)
                    end do
                    B(j,i) = (B(j,i) - innersum)/A(j,j)
                end if
            end do
        end do
    end function ludecomp

end module numericalmethods
