module c2m_basis
use types
implicit none

private
public :: partial_waves_2_operators

integer, parameter :: n_st_oper = 10
! real :: 
contains

subroutine partial_waves_2_operators(pw_param,oper_param)
    implicit none
    real(dp), intent(in)  :: pw_param(:,:)   !< array containing the potential parameters in partial wave basis
    real(dp), intent(out) :: oper_param(:,:) !< array containing the potential parameters in operator basis

    real(dp), allocatable :: spin_isospin_param(:,:)

    call partial_waves_2_spin_isospin(pw_param,spin_isospin_param)
    !call spin_isospin_2_operators(spin_isospin_param,oper_param)
    
end subroutine partial_waves_2_operators

subroutine partial_waves_2_spin_isospin(pw_param,st_param)
    implicit none
    real(dp), intent(in) :: pw_param(:,:) !< array containing the potential parameters in partial wave basis
    real(dp), intent(out), allocatable :: st_param(:,:) !< array containing the potential parameters in spin-isospin basis

    integer :: pw_param_shape(1:2), n_lambdas, n_waves
    
    pw_param_shape = shape(pw_param)
    n_lambdas = pw_param_shape(1)
    n_waves = pw_param_shape(2)

    if(allocated(st_param)) deallocate(st_param)
    allocate(st_param(1:n_lambdas,1:n_st_oper))

end subroutine partial_waves_2_spin_isospin
    
end module c2m_basis