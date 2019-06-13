module c2m_transform
use types, only: dp
use special, only: spherical_bessel_jn
use constants, only: pi
use av18, only: av18_oper_basis, n_operators
use pion_exchange, only: v_one_pion_exch, v_two_pion_exch_nlo, v_two_pion_exch_n2lo
implicit none
private
public delta_shell_2_momentum,sample_av18,transform_all_oper,n_q_operators,n_operators,sample_pion_tail

integer, parameter :: n_q_operators = 24

contains

subroutine transform_all_oper(momentum,lambdas,radii,V_momentum)
    implicit none
    real(dp), intent(in) :: momentum,lambdas(:,:),radii(:)
    real(dp), intent(out) :: V_momentum(:)

    integer :: lambdas_shape(1:2)

    if (size(V_momentum).ne.n_q_operators) then
        print*, 'V_momentum is not the correct size in transform_all_oper'
        stop
    endif

    lambdas_shape = shape(lambdas)

    if (lambdas_shape(2).ne.18) then
         print*, 'lambdas needs to have 18 operators in transform_all_oper'
         stop
     endif 

    if (lambdas_shape(1).ne.size(radii)) then
        print*, 'radii and lambdas(:,i) need to have same size in transform_all_oper'
        stop
    endif

    V_momentum( 1) = delta_shell_2_momentum(momentum,lambdas(:, 1),radii,0,2)
    V_momentum( 2) = delta_shell_2_momentum(momentum,lambdas(:, 2),radii,0,2)
    V_momentum( 3) = delta_shell_2_momentum(momentum,lambdas(:, 3),radii,0,2)
    V_momentum( 4) = delta_shell_2_momentum(momentum,lambdas(:, 4),radii,0,2)
    V_momentum( 5) = delta_shell_2_momentum(momentum,lambdas(:, 5),radii,2,2)
    V_momentum( 6) = delta_shell_2_momentum(momentum,lambdas(:, 6),radii,2,2)
    V_momentum( 7) = delta_shell_2_momentum(momentum,lambdas(:, 7),radii,1,3)
    V_momentum( 8) = delta_shell_2_momentum(momentum,lambdas(:, 8),radii,1,3)
    V_momentum( 9) = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,2,4)
    V_momentum(10) = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,1,3)
    V_momentum(11) = delta_shell_2_momentum(momentum,lambdas(:,10),radii,2,4)
    V_momentum(12) = delta_shell_2_momentum(momentum,lambdas(:,10),radii,1,3)
    V_momentum(13) = delta_shell_2_momentum(momentum,lambdas(:,11),radii,2,4)
    V_momentum(14) = delta_shell_2_momentum(momentum,lambdas(:,11),radii,1,3)
    V_momentum(15) = delta_shell_2_momentum(momentum,lambdas(:,12),radii,2,4)
    V_momentum(16) = delta_shell_2_momentum(momentum,lambdas(:,12),radii,1,3)
    V_momentum(17) = delta_shell_2_momentum(momentum,lambdas(:,13),radii,2,4)
    V_momentum(18) = delta_shell_2_momentum(momentum,lambdas(:,13),radii,1,3)
    V_momentum(19) = delta_shell_2_momentum(momentum,lambdas(:,14),radii,2,4)
    V_momentum(20) = delta_shell_2_momentum(momentum,lambdas(:,14),radii,1,3)
    V_momentum(21) = delta_shell_2_momentum(momentum,lambdas(:,15),radii,0,2)
    V_momentum(22) = delta_shell_2_momentum(momentum,lambdas(:,16),radii,0,2)
    V_momentum(23) = delta_shell_2_momentum(momentum,lambdas(:,17),radii,2,2)
    V_momentum(24) = delta_shell_2_momentum(momentum,lambdas(:,18),radii,0,2)

    V_momentum = V_momentum!/(2*pi)**3
    
end subroutine transform_all_oper

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

subroutine sample_pion_tail(radii,lecs,tail)
    implicit none
    real(dp), intent(in) :: radii(:),lecs(:)
    real(dp), intent(out) :: tail(:,:)

    integer :: tail_shape(1:2),i,n_lambdas
    real(dp), allocatable :: temp_array(:,:)
    real(dp) :: v_ope(1:n_operators), v_tpe_nlo(1:n_operators), v_tpe_n2lo(1:n_operators)
    real(dp) :: r_i, delta_r

    tail_shape = shape(tail)
    allocate(temp_array(1:tail_shape(2),1:tail_shape(1)))

    if(size(radii).ne.tail_shape(1)) then
        print*, 'incompatible sizes between radii and tail in sample_pion_tail'
        stop
    endif

    if (tail_shape(2).ne.n_operators) then
        print*, 'incorrect size for tail in sample_pion_tail'
        stop
    endif

    n_lambdas = tail_shape(1)
    do i=1,n_lambdas
        r_i = radii(i)
        call v_one_pion_exch(r_i,v_ope)
        call v_two_pion_exch_nlo(r_i,v_tpe_nlo)
        call v_two_pion_exch_n2lo(r_i,lecs,v_tpe_n2lo)
        temp_array(:,i) = v_ope + v_tpe_nlo + v_tpe_n2lo
    enddo
    delta_r = radii(2) - radii(1)
    temp_array = temp_array*delta_r
    tail = transpose(temp_array)

end subroutine sample_pion_tail

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

    r = 4*r*pi/q**q_power
end function delta_shell_2_momentum
    
end module c2m_transform