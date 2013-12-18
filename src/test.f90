program test
    use conductionfem
    implicit none
	real(8) :: t1,t2
	call cpu_time(t1)
    call finitesolution()
	call cpu_time(t2)
	write(*,*) "Total time: ", t2-t1
end program test
