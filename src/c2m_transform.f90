module c2m_transform
use types, only: dp
use special, only: spherical_bessel_jn
use constants, only: pi
use av18, only: av18_oper_basis, n_operators
implicit none
private
public delta_shell_2_momentum

contains

subroutine sample_av18(delta_r,r_max,radii,av18_lambdas)
    implicit none
    real(dp), intent(in) :: delta_r,r_max
    real(dp), intent(out), allocatable :: radii(:),av18_lambdas(:,:)
    integer :: n_lambdas, i
    real(dp) :: v_oper(1:n_operators), r_i
    real(dp), allocatable :: temp_array(:,:)

    n_lambdas = int(r_max/delta_r)
    if(allocated(radii)) deallocate(radii) 
    if(allocated(av18_lambdas)) deallocate(av18_lambdas)

    allocate(radii(1:n_lambdas))
    allocate(av18_lambdas(1:n_lambdas,1:n_operators))
    allocate(temp_array(1:n_operators,n_lambdas))

    radii = [((i-0.5_dp)*delta_r,i=1,n_lambdas)]

    do i=1,n_lambdas
        r_i = radii(i)
        call av18_oper_basis(r_i,v_oper)
        temp_array(:,i) = v_oper
    enddo
    temp_array = temp_array*delta_r
    av18_lambdas = transpose(temp_array)
end subroutine sample_av18

real(dp) function delta_shell_2_momentum(q,lambdas,radii,q_power,r_power) result(r)
    implicit none
    real(dp), intent(in) :: q, lambdas(:), radii(:)
    integer, intent(in) :: q_power, r_power

    integer :: n_lambdas,i
    real(dp) :: r_i, lambda_i

    n_lambdas = size(lambdas)
    if (n_lambdas.ne.size(radii)) then
        print*, 'lambdas and radii have incompatible siez in delta_shell_2_momentum'
        stop    
    endif

    r = 0._dp
    do i = 1, n_lambdas
        r_i = radii(i)
        lambda_i = lambdas(i)
        r = r + spherical_bessel_jn(q_power,q*r_i)*lambda_i*r_i**r_power
    enddo

    r = r*pi/q**q_power
end function delta_shell_2_momentum
    
end module c2m_transform