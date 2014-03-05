program test
	use htfem
    implicit none
	integer :: values1(8),values2(8)
	character(8)  :: date
	character(10) :: time
	character(5)  :: zone
	call date_and_time(date,time,zone,values1)
	write(*,'(a,(3i8,2x),i8)') "Start time: ", values1(5:8)
	call fem()
	call date_and_time(date,time,zone,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)
	write(*,'(a,(3i8,2x),i8)') "Runtime: ",values2(5:8)-values1(5:8)
end program test
