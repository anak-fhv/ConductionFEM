! some auxialiary functions for raytracing
! does not require any other module

module rt_helper

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
    
end module rt_helper