module c2m_montecarlo
use types, only: dp
use statistics, only: mean, sample_variance
use c2m_basis, only: allocate_operators, partial_waves_2_operators
use c2m_transform, only: delta_shell_2_momentum, transform_all_oper
implicit none

private
public mc_momentum_dependence

contains

subroutine mc_momentum_dependence(param_samples,momentum,dr,V_mean,V_variance)
    implicit none
    real(dp), intent(in) :: param_samples(:,:,:)
    real(dp), intent(in) :: momentum, dr
    real(dp), intent(out) :: V_mean(:), V_variance(:)

    integer, parameter :: n_operators=24

    integer :: samples_shape(1:3), n_samples, n_lambdas, n_waves
    real(dp), allocatable :: oper_parameters(:,:), pw_parameters(:,:), radii(:)
    real(dp), allocatable :: temp_array(:,:), V_samples(:,:), V_momentum(:)
    integer :: i, n_q_operators

    samples_shape = shape(param_samples)
    n_lambdas = samples_shape(1)
    n_samples = samples_shape(3)
    n_q_operators = size(V_mean)

    allocate(pw_parameters,mold=param_samples(:,:,1))
    allocate(radii(1:n_lambdas))
    allocate(V_momentum(1:n_q_operators))
    allocate(temp_array(1:n_q_operators,1:n_samples))
    allocate(V_samples(1:n_samples,1:n_q_operators))
    call allocate_operators(n_lambdas,oper_parameters)

    radii = [(i*dr,i=1,n_lambdas)]

    do i=1,n_samples
        oper_parameters = 0
        pw_parameters = param_samples(:,:,i)
        call partial_waves_2_operators(pw_parameters,oper_parameters)
        call transform_all_oper(momentum,oper_parameters(:,1:18),radii,V_momentum)
        temp_array(:,i) = V_momentum 
    enddo
    V_samples = transpose(temp_array)
    do i = 1,n_q_operators
        V_mean(i) = mean(V_samples(:,i))
        V_variance(i) = sample_variance(V_samples(:,i))
    enddo
    
end subroutine mc_momentum_dependence
    
end module c2m_montecarlo