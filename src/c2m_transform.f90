module c2m_transform
use types, only: dp
use special, only: spherical_bessel_jn
use constants, only: pi
use av18, only: av18_oper_basis, n_operators
use pion_exchange, only: v_one_pion_exch, v_two_pion_exch_nlo, v_two_pion_exch_n2lo
implicit none
private
public delta_shell_2_momentum,sample_av18,transform_all_oper,n_q_operators,n_operators,sample_pion_tail,n_av18q_operators

integer, parameter :: n_q_operators = 12
integer, parameter :: n_av18q_operators = 24

contains

subroutine transform_all_oper(momentum,lambdas,radii,V18q,V_momentum)
    implicit none
    real(dp), intent(in) :: momentum,lambdas(:,:),radii(:)
    real(dp), intent(out) :: V18q(:), V_momentum(:)

    integer :: lambdas_shape(1:2)
    real(dp) :: p_momentum, q_dot_p

    if (size(V_momentum).ne.n_q_operators) then
        print*, 'V_momentum is not the correct size in transform_all_oper'
        stop
    endif

    if (size(V18q).ne.n_av18q_operators) then
        print*, 'V18q is not the correct size in transform_all_oper'
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

    V18q( 1) = delta_shell_2_momentum(momentum,lambdas(:, 1),radii,0,2)
    V18q( 2) = delta_shell_2_momentum(momentum,lambdas(:, 2),radii,0,2)
    V18q( 3) = delta_shell_2_momentum(momentum,lambdas(:, 3),radii,0,2)
    V18q( 4) = delta_shell_2_momentum(momentum,lambdas(:, 4),radii,0,2)
    V18q( 5) = delta_shell_2_momentum(momentum,lambdas(:, 5),radii,2,2)
    V18q( 6) = delta_shell_2_momentum(momentum,lambdas(:, 6),radii,2,2)
    V18q( 7) = delta_shell_2_momentum(momentum,lambdas(:, 7),radii,1,3)
    V18q( 8) = delta_shell_2_momentum(momentum,lambdas(:, 8),radii,1,3)
    V18q( 9) = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,2,4)
    V18q(10) = delta_shell_2_momentum(momentum,lambdas(:, 9),radii,1,3)
    V18q(11) = delta_shell_2_momentum(momentum,lambdas(:,10),radii,2,4)
    V18q(12) = delta_shell_2_momentum(momentum,lambdas(:,10),radii,1,3)
    V18q(13) = delta_shell_2_momentum(momentum,lambdas(:,11),radii,2,4)
    V18q(14) = delta_shell_2_momentum(momentum,lambdas(:,11),radii,1,3)
    V18q(15) = delta_shell_2_momentum(momentum,lambdas(:,12),radii,2,4)
    V18q(16) = delta_shell_2_momentum(momentum,lambdas(:,12),radii,1,3)
    V18q(17) = delta_shell_2_momentum(momentum,lambdas(:,13),radii,2,4)
    V18q(18) = delta_shell_2_momentum(momentum,lambdas(:,13),radii,1,3)
    V18q(19) = delta_shell_2_momentum(momentum,lambdas(:,14),radii,2,4)
    V18q(20) = delta_shell_2_momentum(momentum,lambdas(:,14),radii,1,3)
    V18q(21) = delta_shell_2_momentum(momentum,lambdas(:,15),radii,0,2)
    V18q(22) = delta_shell_2_momentum(momentum,lambdas(:,16),radii,0,2)
    V18q(23) = delta_shell_2_momentum(momentum,lambdas(:,17),radii,2,2)
    V18q(24) = delta_shell_2_momentum(momentum,lambdas(:,18),radii,0,2)

    call av18_momentum_2_spin_isospin(p_momentum,momentum,q_dot_p,V18q,V_momentum)
    
end subroutine transform_all_oper

subroutine av18_momentum_2_spin_isospin(p,q,q_dot_p,V18q,Vstq)
    implicit none
    real(dp), intent(in) :: p, q, q_dot_p, V18q(:)
    real(dp), intent(out) :: Vstq(:)

    real(dp) :: V_1, V_2, V_3, V_4, V_5, V_6, V_7, V_8, V_9a, V_9b, V_10a, V_10b, V_11a, &
        & V_11b, V_12a, V_12b, V_13a, V_13b, V_14a, V_14b, V_15, V_16, V_17, V_18

    if (size(Vstq).ne.n_q_operators) then
        print*, 'Vstq is not the correct size in av18_momentum_2_spin_isospin'
        stop
    endif

    if (size(V18q).ne.n_av18q_operators) then
        print*, 'V18q is not the correct size in av18_momentum_2_spin_isospin'
        stop
    endif

    V_1   = V18q( 1)
    V_2   = V18q( 2)
    V_3   = V18q( 3)
    V_4   = V18q( 4)
    V_5   = V18q( 5)
    V_6   = V18q( 6)
    V_7   = V18q( 7)
    V_8   = V18q( 8)
    V_9a  = V18q( 9)
    V_9b  = V18q(10)
    V_10a = V18q(11)
    V_10b = V18q(12)
    V_11a = V18q(13)
    V_11b = V18q(14)
    V_12a = V18q(15)
    V_12b = V18q(16)
    V_13a = V18q(17)
    V_13b = V18q(18)
    V_14a = V18q(19)
    V_14b = V18q(20)
    V_15  = V18q(21)
    V_16  = V18q(22)
    V_17  = V18q(23)
    V_18  = V18q(24)


    Vstq( 1) = V_1 - (V_9a + V_10a + V_11a + V_12a + V_13a/2)*&
        & (p**2*q**2 - q_dot_p**2) + (2+V_9b + V_13b)*(p**2 - q**2/4)
    Vstq( 2) = V_2 - (V_12a + V_14a/2)*(p**2*q**2 - q_dot_p**2) &
        & + (2*V_10b + V_14b)*(p**2 - q**2/4)
    Vstq( 3) = V_3 + (2*V_11b + V_13b)*p**2 - (V_11b/2 + V_13b/4 - V_5)*q**2
    Vstq( 4) = V_4 + (2*V_12b + V_14b)*(p**2 - q**2/4) + V_6*q**2
    Vstq( 5) = V_13b/8 - 3*V_5
    Vstq( 6) = V_14b/8 - 3*V_6
    Vstq( 7) = V_7
    Vstq( 8) = V_8
    Vstq( 9) = - V_13a/2
    Vstq(10) = - V_14a/2
    Vstq(11) = - V_13b/2
    Vstq(12) = - V_14b/2

    
end subroutine av18_momentum_2_spin_isospin

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

subroutine sample_pion_tail(radii,lecs,chiral_order,tail)
    implicit none
    real(dp), intent(in) :: radii(:),lecs(:)
    integer, intent(in) :: chiral_order
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
        if(chiral_order .ge. 0) then
            call v_one_pion_exch(r_i,v_ope)
        else
            v_ope = 0._dp
        endif
        if (chiral_order .ge. 1) then
            call v_two_pion_exch_nlo(r_i,v_tpe_nlo)
        else
            v_tpe_nlo = 0._dp
        endif
        if (chiral_order .ge. 2) then
            call v_two_pion_exch_n2lo(r_i,lecs,v_tpe_n2lo)
        else
            v_tpe_n2lo = 0._dp
        endif
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