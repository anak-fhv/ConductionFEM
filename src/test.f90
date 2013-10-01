program test
    use readmesh
    use numericalmethods
    use finitesolver
    implicit none
!    double precision, dimension(:,:), allocatable :: vertices, stiffness, globalstiffness
!    integer, dimension(:), allocatable :: surfacenodes, domainelements
!    double precision, dimension(:), allocatable :: elementvolumes
!    character, dimension(:), allocatable :: surfacenames
!    integer, dimension(7) :: meshdetails
    integer, dimension(:,:), allocatable :: connectivity, surfacefaces
    double precision, dimension(3,3) :: A,B,C
!    double precision, dimension(4,4) :: A
!    double precision, dimension(2) :: timearray
!    double precision :: dm, result, T1, T2
    integer :: i
    call finitesolution()
end program test
