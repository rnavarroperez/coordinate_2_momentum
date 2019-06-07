module c2m_montecarlo
use types, only: dp
use statistics, only: mean, sample_variance
use c2m_basis, only: allocate_operators, partial_waves_2_operators
use c2m_transform, only: delta_shell_2_momentum
implicit none

private
public mc_momentum_dependence

contains

subroutine mc_momentum_dependence(param_samples,dr)
    implicit none
    real(dp), intent(in) :: param_samples(:,:,:)
    real(dp), intent(in) :: dr

    integer :: samples_shape(1:3), n_samples, n_lambdas, n_waves
    real(dp), allocatable :: oper_parameters(:,:), pw_parameters(:,:), radii(:),lambdas14(:)
    real(dp), allocatable :: Vmomentum14a(:), Vmomentum14b(:)
    integer :: i
    real(dp) :: momentum

    samples_shape = shape(param_samples)
    n_lambdas = samples_shape(1)
    n_samples = samples_shape(3)

    allocate(pw_parameters,mold=param_samples(:,:,1))
    allocate(radii(1:n_lambdas))
    allocate(lambdas14,mold=radii)
    allocate(Vmomentum14a(1:n_samples))
    allocate(Vmomentum14b,mold=Vmomentum14a)
    call allocate_operators(n_lambdas,oper_parameters)

    radii = [(i*dr,i=1,n_lambdas)]
    momentum = 1.e-1_dp
    
    do i=1,n_samples
        oper_parameters = 0
        pw_parameters = param_samples(:,:,i)
        call partial_waves_2_operators(pw_parameters,oper_parameters)
        lambdas14 = oper_parameters(:,14)
        Vmomentum14a(i) = delta_shell_2_momentum(momentum,lambdas14,radii,2,4)
        Vmomentum14b(i) = delta_shell_2_momentum(momentum,lambdas14,radii,1,3)
    enddo
    print*, mean(Vmomentum14a), sample_variance(Vmomentum14a)
    print*, mean(Vmomentum14b), sample_variance(Vmomentum14b)

    
end subroutine mc_momentum_dependence
    
end module c2m_montecarlo