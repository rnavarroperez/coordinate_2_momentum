module c2m_montecarlo
use types, only: dp
use statistics, only: mean, sample_variance
use c2m_basis, only: allocate_operators, partial_waves_2_operators
use c2m_transform, only: delta_shell_2_momentum, transform_all_oper,n_av18_operators=>n_operators,sample_pion_tail
implicit none

private
public mc_momentum_dependence

contains

subroutine mc_momentum_dependence(param_samples,momentum,dr,dr_tail,r_max,V_mean,V_variance)
    implicit none
    real(dp), intent(in) :: param_samples(:,:,:)
    real(dp), intent(in) :: momentum, dr,dr_tail,r_max
    real(dp), intent(out) :: V_mean(:), V_variance(:)

    integer, parameter :: n_lecs=3

    integer :: samples_shape(1:3), n_samples, n_lambdas, n_waves,n_tail
    real(dp), allocatable :: oper_parameters(:,:), pw_parameters(:,:), radii(:),r_tail(:)
    real(dp), allocatable :: temp_array(:,:), V_samples(:,:), V_momentum(:)
    real(dp), allocatable :: tail_oper_basis(:,:)
    real(dp) :: rcut, lecs(1:n_lecs)
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
    rcut = n_lambdas*dr
    n_tail = int(ceiling((r_max-rcut)/dr_tail))
    r_tail = [(rcut + (i-0.5_dp)*dr_tail,i=1,n_tail)]

    allocate(tail_oper_basis(1:n_tail,1:n_av18_operators))
    
    do i=1,n_samples
        oper_parameters = 0
        pw_parameters = param_samples(:,:,i)
        lecs = pw_parameters(1:n_lecs,19)
        call partial_waves_2_operators(pw_parameters,oper_parameters)
        call transform_all_oper(momentum,oper_parameters(:,1:n_av18_operators),radii,V_momentum)
        temp_array(:,i) = V_momentum 
        !call sample_pion_tail(r_tail,lecs,tail_oper_basis)
        !call transform_all_oper(momentum,tail_oper_basis,r_tail,V_momentum)
        !temp_array(:,i) = temp_array(:,i) + V_momentum 
    enddo
    V_samples = transpose(temp_array)
    do i = 1,n_q_operators
        V_mean(i) = mean(V_samples(:,i))
        V_variance(i) = sample_variance(V_samples(:,i))
    enddo
    
end subroutine mc_momentum_dependence
    
end module c2m_montecarlo