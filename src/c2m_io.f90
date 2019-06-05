module c2m_io

use types

implicit none

private
public :: read_parameters

character(len=25), parameter :: f_name = 'tpe218pwanew.dat' !< file nane subfix

contains

subroutine read_parameters(param)
    implicit none
    real(dp), intent(out) :: param(:,:) !< array to be filled with the potential parameters in partial wave basis

    integer            :: unit, ierr,readerr,i_waves,i_lambdas
    character(len=128) :: filename
    logical            :: exists
    integer            :: param_shape(1:2),n_lambdas,n_waves

    param_shape = shape(param)
    n_lambdas = param_shape(1)
    n_waves = param_shape(2)

    filename = 'BestLambdas'//f_name
    inquire(file=trim(filename), exist=exists)
    if (exists) then
        open(newunit = unit, file=trim(filename), action="read", iostat=ierr)
        if (ierr .eq. 0) then
            do i_waves=1,n_waves
                read(unit, *, iostat=readerr) param(:,i_waves)
                if (readerr .ne. 0) then
                    print*, "Error ", readerr, & 
                            " attempting to read file ", & 
                            trim(filename), "In wave number", i_waves
                    stop
                endif
                read(unit, *, iostat=readerr) 
            enddo
        else
            print*, "Error ", ierr ," attempting to open file ", trim(filename)
            stop
        endif
    else
        print*, "Error -- cannot find file: ", trim(filename)
        stop
    endif

end subroutine read_parameters

    
end module c2m_io