module c2m_io

use types
use c2m_transform
implicit none

private
public :: read_parameters, write_momentum_dependence

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
        close(unit)
    else
        print*, "Error -- cannot find file: ", trim(filename)
        stop
    endif

end subroutine read_parameters

subroutine write_momentum_dependence(oper_parameters,dr)
    implicit none
    real(dp), intent(in) :: oper_parameters(:,:), dr
    integer :: unit,ierr
    character(len=128) :: filename
    real(dp) :: momentum,V_mom_14a,V_mom_14b
    real(dp), allocatable :: lambdas14(:), radii(:)
    integer :: oper_parameters_shape(1:2),n_lambdas,i

    oper_parameters_shape = shape(oper_parameters)
    n_lambdas = oper_parameters_shape(1)

    allocate(lambdas14(1:n_lambdas),radii(1:n_lambdas))

    lambdas14 = oper_parameters(:,14)
    do i = 1,n_lambdas
        radii(i) = i*dr
    enddo

    filename = 'momentum_dependence'//f_name

    open(newunit = unit, file=trim(filename), action="write", iostat=ierr)
    if (ierr .eq. 0) then
        momentum = 1.e-1_dp
        do
            V_mom_14a = delta_shell_2_momentum(momentum,lambdas14,radii,2,4)
            V_mom_14b = delta_shell_2_momentum(momentum,lambdas14,radii,1,3)
            write(unit,'(4e21.8)') momentum, V_mom_14a, V_mom_14b
            momentum = momentum + 1.e-1_dp
            if (momentum.gt.100._dp) exit
        enddo
        
    else
        print*, "Error ", ierr ," attempting to open file ", trim(filename)
        stop
    endif
    close(unit)
    
end subroutine write_momentum_dependence
    
end module c2m_io