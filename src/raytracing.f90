! RayTracing
! author: Steffen Finck
! contact: steffen.finck@fhv.at
! date: 21.02.2014
! version: 0.1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

program raytracing

    use math_funs
    use pre_process_data
    implicit none
    
    real(dp), dimension(3)             :: t1, t2, t3, t4 ! points of the tetraeder
    real(dp), dimension(3,4)           :: T, vp
    real(dp), dimension(3)             :: ro, rd  ! origin and direction of ray
    integer, dimension(12), parameter :: sorder = [1, 2, 3, 2, 3, 4, 1, 4, 3, 1, 2, 4] 
    integer                           :: i, n = 1000000
    integer                           :: nvp, c1, c2, nvp_total, c1_total, c2_total ! vectors and points 
                                                                                 ! for each surface 
    real                              :: time1, time2 ! runtime variables
    character(len=20)                          :: file_name
    type(tetraElement), dimension(:), allocatable :: tetras, tetras1
    integer                         :: npart
    
    ! user input
    file_name = "ex10.msh"
    npart = 4
    call read_mesh_data(file_name, npart, tetras)
    
!     write(*,*) tetras(1)
!     write(*,*)
    
!     npart = 1
!     call read_mesh_data(file_name, npart, tetras1)
!     write(*,*) tetras(1)
    !write(*,*) vert
 
!  ! initialize parameter
!  t1 = [-1d0, 1d0, -1d0]
!  t2 = [0d0 , -1d0, -1d0]
!  t3 = [1d0 , 0d0, -1d0]
!  t4 = [0d0 , 0d0, 1d0]
!  T = reshape((/t1,t2,t3,t4/),(/3,4/))
!  
!  !write(*,*) t1
!  !write(*,*) t2
!  !write(*,*) t3
!  !write(*,*) t4
!  !write(*,*) T
!  !write(*,*) T(2,2)
! 
!  c1_total = 0
!  c2_total = 0
!  nvp_total = 0
!  
!  call cpu_time(time1)
!  
!  do i = 1,n
!  
!   ! create ray
!   ro = [2d0*drand(0), 2d0*drand(0), 2d0*drand(0)]
!   rd = -2d0*ro
!   rd = rd/norm(rd)
!  
!   !write(*,*)
!   !write(*,*) ro
!   !write(*,*) rd
!  
!   call singleTetra(T, ro, rd, sorder, vp, nvp, c1, c2)
!   c1_total = c1_total + c1
!   c2_total = c2 + c2_total
!   if ((nvp /= 2) .and. (nvp /= 0)) then
!    nvp_total = nvp_total + 1  
!   end if 
!  end do
!  
!  call cpu_time(time2)
!  
!  write(*,*) "Time: ", time2 - time1
!  write(*,*) "Total number of skipped 1st cycles: ", c1_total
!  write(*,*) "Total number of skipped 2nd cycles: ", c2_total
!  write(*,*) "Relative number of wrong results: ", real(nvp_total)/real(n) 
 
end program raytracing 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! SUBROUTINES 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! ! subroutine to find intersection point for single tetraeder 
! subroutine singleTetra(T, ro, rd, sorder, vp, nvp, c1, c2)
!  
!  use mathfuns
!  implicit none
!  
!  real(8), dimension(3,4), intent(in)  :: T           ! matrix with tetraeder points
!  real(8), dimension(3), intent(in)    :: ro, rd      ! origin and direction of ray
!  integer, dimension(12), intent(in)   :: sorder      ! array containing the 3 points for each 
!                                                      ! of the 4 surfaces of the tetraeder
!  real(8), dimension(3,2), intent(out) :: vp          ! matrix of valid points
!  real(8), dimension(3,4)              :: vp_temp
!  integer, intent(out)                 :: nvp, c1, c2 ! some stats about the routine
!  real(8), dimension(3)             :: cg             ! center of gravity
!  real(8)                           :: m_dist         ! maximal distance of the tetraeder points 
!                                                      ! from the center of gravity
!  integer                           :: s, id            ! counter variables
!  real(8), dimension(3)             :: rp, v1, v2, nsf, p ! vectors and points for each surface 
!  real(8)                           :: d, alpha, dp       ! surface parameters
!  integer, dimension(4)             :: vp_id 
!  real(8), dimension(3)             :: rhs
!  real(8)                           :: beta1, beta2
!  
!  ! get center of gravity
!  cg = 1/4d0*(T(:,1) + T(:,2) + T(:,3) + T(:,4))
!  
!  !write(*,*)
!  !write(*,*) cg
! 
!  ! get distance of the tetraeder points to the enter of gravity 
!  m_dist = maxval([norm(T(:,1)-cg), norm(T(:,2)-cg), norm(T(:,3)-cg), norm(T(:,4)-cg)])
!  
!  ! loop over surfaces
!  nvp = 0
!  c1 = 4
!  c2 = 4
!  vp_temp = -10
!  vp_id = 0
!  do s = 1,4
!   
!   rp = T(:,sorder(3*s - 2))       ! reference point for surface
!   v1 = T(:,sorder(3*s - 1)) - rp  ! first vector defining surface
!   v2 = T(:,sorder(3*s)) - rp      ! 2nd vector defining surface
!   
!   ! get normal vector of surface 
!   nsf = cross(v1,v2)
!   nsf = nsf/norm(nsf) ! normalize
!   
!   ! ray direction and normal vector are parallel
!   dp = dot_product(nsf, rd)
!   if (abs(dp) <= 1d-14) cycle
!   c1 = c1 - 1
!   
!   ! calculate intersection point p between ray and surface
!   d = dot_product(nsf,rp)  ! surface parameter
!   alpha = (d - dot_product(nsf,ro))/dp
!   p = ro + alpha*rd;
!   
!   if (norm(p - cg) > m_dist) cycle  ! cycle has a distance than largest distance in tetraeder
!   c2 = c2 - 1
!   integer, parameter    :: dp=selected_real_kind(p=14)
!   ! point could be considered as valid intersection point
!   !write(*,*) s
!   !write(*,*) p
!   nvp = nvp + 1;
!   vp_temp(:,s) = p
!   vp_id(s) = 1
!  end do
! 
!  !write(*,*) vp
!  !write(*,*)
!  !write(*,*) vp_id
!   
!  ! if there are not the right number of intersection points
!  ! one has to look for the intersection points on the surface of the
!  ! tetraeder 
!  
!  if (nvp /= 2) then
!   do s = 1,4
!    if (vp_id(s) /= 1) cycle
!    
!    vp_id(s) = 0
!    
!    rp = T(:,sorder(3*s - 2))       ! reference point for surface
!    v1 = T(:,sorder(3*s - 1)) - rp  ! first vector defining surface
!    v2 = T(:,sorder(3*s)) - rp      ! 2nd vector defining surface
!    
!    p = vp_temp(:,s)
!    
! ! !    if (norm(p-v1+rp) <1d-14) then
! ! !     vp_id(s) = 1
! ! !     cycle
! ! !    else if (norm(p-v2+rp) <1d-14) then
! ! !     vp_id(s) = 1
! ! !     cycle 
! ! !    end if
!    
!    rhs = p - rp
!    
!    beta2 = (v1(2)*rhs(1) - rhs(2)*v1(1))/(v1(2)*v2(1) - v2(2)*v1(1))
!    beta1 = (rhs(1) - v2(1)*beta2)/v1(1)
!    
!    !!!! above still needs checking that no division by zero happens!!!!
!    
!    if (abs(v1(3)*beta1 + v2(3)*beta2 - rhs(3)) > 1d-14) then
!     write(*,*) "something went wrong!"
!     write(*,*) rhs(3) - (v1(3)*beta1 + v2(3)*beta2)
!     return
!    end if
!     
!    if ((beta1 + beta2 <= 1) .and. (beta1 >= 0) .and. (beta2 >= 0)) vp_id(s) = 1  
!    
!    if (sum(vp_id) == 2) exit
! 
!   end do
!  end if 
!  
!  id = 1;
!  do s=1,4
!   if (vp_id(s) == 1) then
!    vp(:,id) = vp_temp(:,s)
!    id = id + 1
!   end if
!  end do
!  
!  nvp = id-1
!  
!  return
!  
! end subroutine singleTetra