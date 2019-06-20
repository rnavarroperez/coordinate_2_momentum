module c2m_transform
use types, only: dp
use special, only: spherical_bessel_jn
use constants, only: pi
use av18, only: av18_oper_basis, n_operators
use pion_exchange, only: v_one_pion_exch, v_two_pion_exch_nlo, v_two_pion_exch_n2lo
implicit none
private
public delta_shell_2_momentum,sample_av18,transform_all_oper,n_q_operators,n_operators,sample_pion_tail

integer, parameter :: n_q_operators = 12

contains

subroutine transform_all_oper(momentum,lambdas,radii,V_momentum)
    implicit none
    real(dp), intent(in) :: momentum,lambdas(:,:),radii(:)
    real(dp), intent(out) :: V_momentum(:)

    integer :: lambdas_shape(1:2)
    real(dp) :: V_1, V_2, V_3, V_4, V_5, V_6, V_7, V_8, V_9a, V_9b, V_10a, V_10b, V_11a, &
        & V_11b, V_12a, V_12b, V_13a, V_13b, V_14a, V_14b, V_15, V_16, V_17, V_18, &
        & p_momentum, q_dot_p

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

    p_momentum = 0._dp
    q_dot_p = 0._dp

    V_1   = delta_shell_2_momentum(momentum,lambdas(:, 1),radii,0,2)
    V_2   = delta_shell_2_momentum(momentum,lambdas(:, 2),radii,0,2)
    V_3   = delta_shell_2_momentum(momentum,lambdas(:, 3),radii,0,2)
    V_4   = delta_shell_2_momentum(momentum,lambdas(:, 4),radii,0,2)
    V_5   = delta_shell_2_momentum(momentum,lambdas(:, 5),radii,2,2)
    V_6   = delta_shell_2_momentum(momentum,lambdas(:, 6),radii,2,2)
    V_7   = delta_shell_2_momentum(momentum,lambdas(:, 7),radii,1,3)
    V_8   = delta_shell_2_momentum(momentum,lambdas(:, 8),radii,1,3)
    V_9a  = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,2,4)
    V_9b  = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,1,3)
    V_10a = delta_shell_2_momentum(momentum,lambdas(:,10),radii,2,4)
    V_10b = delta_shell_2_momentum(momentum,lambdas(:,10),radii,1,3)
    V_11a = delta_shell_2_momentum(momentum,lambdas(:,11),radii,2,4)
    V_11b = delta_shell_2_momentum(momentum,lambdas(:,11),radii,1,3)
    V_12a = delta_shell_2_momentum(momentum,lambdas(:,12),radii,2,4)
    V_12b = delta_shell_2_momentum(momentum,lambdas(:,12),radii,1,3)
    V_13a = delta_shell_2_momentum(momentum,lambdas(:,13),radii,2,4)
    V_13b = delta_shell_2_momentum(momentum,lambdas(:,13),radii,1,3)
    V_14a = delta_shell_2_momentum(momentum,lambdas(:,14),radii,2,4)
    V_14b = delta_shell_2_momentum(momentum,lambdas(:,14),radii,1,3)
    V_15  = delta_shell_2_momentum(momentum,lambdas(:,15),radii,0,2)
    V_16  = delta_shell_2_momentum(momentum,lambdas(:,16),radii,0,2)
    V_17  = delta_shell_2_momentum(momentum,lambdas(:,17),radii,2,2)
    V_18  = delta_shell_2_momentum(momentum,lambdas(:,18),radii,0,2)

    V_momentum( 1) = V_1 - (V_9a + V_10a + V_11a + V_12a + V_13a/2)*&
        & (p_momentum**2*momentum**2 - q_dot_p**2) + (2+V_9b + V_13b)*&
        & (p_momentum**2 - momentum**2/4)
    V_momentum( 2) = V_2 - (V_12a + V_14a/2)*(p_momentum**2*momentum**2 - q_dot_p**2) &
        & + (2*V_10b + V_14b)*(p_momentum**2 - momentum**2/4)
    V_momentum( 3) = V_3 + (2*V_11b + V_13b)*p_momentum**3 &
        & - (V_11b/2 + V_13b/4 - V_5)*momentum**2
    V_momentum( 4) = V_4 + (2*V_12b + V_14b)*(p_momentum**2 - momentum**2/4) &
        & + V_6*momentum**2
    V_momentum( 5) = V_13b/8 - 3*V_5
    V_momentum( 6) = V_14b/8 - 3*V_6
    V_momentum( 7) = V_7
    V_momentum( 8) = V_8
    V_momentum( 9) = - V_13a/2
    V_momentum(10) = - V_14a/2
    V_momentum(11) = - V_13b/2
    V_momentum(12) = - V_14b/2
    
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