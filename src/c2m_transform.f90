module c2m_transform
use types
use special
use constants, only: pi
implicit none
private
public delta_shell_2_momentum

contains

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