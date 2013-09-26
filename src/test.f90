program test
    use readmesh
    use numericalmethods
    use finitesolver
    implicit none
!    real, dimension(:,:), allocatable :: vertices, stiffness, globalstiffness
!    integer, dimension(:), allocatable :: surfacenodes, domainelements
!    real, dimension(:), allocatable :: elementvolumes
!    character, dimension(:), allocatable :: surfacenames
!    integer, dimension(7) :: meshdetails
    integer, dimension(:,:), allocatable :: connectivity, surfacefaces
    real, dimension(3,3) :: A,B
!    real, dimension(4,4) :: A
!    real, dimension(2) :: timearray
!    real :: dm, result, T1, T2
!    integer :: i
    print *, 'testing'
    A = reshape((/10,-3,1,2,-6,1,-1,2,5/),(/3,3/))
    B = ludecomp(A)
    print *, B(1,:)
    print *, B(2,:)
    print *, B(3,:)
    print *, 'testing again'
end program test
