module c2m_basis
use types, only: dp
use constants, only: hbar_c_MeV_fm, proton_mass_MeV
implicit none

private
public :: partial_waves_2_operators, set_fit_flags, allocate_operators

integer, parameter :: n_operators = 21

contains

subroutine allocate_operators(n_lambdas,oper_param)
    implicit none
    integer, intent(in) :: n_lambdas
    real(dp), allocatable, intent(out) :: oper_param(:,:) 

    if(allocated(oper_param)) deallocate(oper_param)
    allocate(oper_param(1:n_lambdas,1:n_operators))
end subroutine allocate_operators

subroutine set_fit_flags(pw_param,flags)
    implicit none
    real(dp), intent(in) :: pw_param(:,:) !< array containing the potential parameters in partial wave basis (in units of fm^{-1})
    logical, intent(out) :: flags(:,:) !< array containing logial flags to indicate which parameters were fitted

    if(any(shape(pw_param).ne.shape(flags))) then
        print*, 'Inconsistent shapes in set_fit_flags'
        stop
    endif

    where (pw_param.eq.0._dp)
        flags = .false.
    else where
        flags = .true.
    end where
    
end subroutine set_fit_flags

subroutine partial_waves_2_operators(pw_param,oper_param)
    implicit none
    real(dp), intent(in)  :: pw_param(:,:)   !< array containing the potential parameters in partial wave basis (in units of fm^{-1})
    real(dp), intent(out) :: oper_param(:,:) !< array containing the potential parameters in operator basis (return in units of MeV*fm)

    real(dp), allocatable, dimension(:) :: l_1s0np, l_1s0pp, l_3p0,    &
    & l_1p1, l_3p1, l_3s1, l_ep1, l_3d1, l_1d2, l_3d2, l_3p2, l_ep2,   &
    & l_3f2, l_1f3, l_3d3
    
    integer :: pw_param_shape(1:2), n_lambdas,i

    pw_param_shape = shape(pw_param)
    n_lambdas = pw_param_shape(1)

    allocate(l_1s0np(1:n_lambdas))
    allocate(l_1s0pp,l_3p0,l_1p1,l_3p1,l_3s1,l_ep1,l_3d1,l_1d2,l_3d2,  &
    &        l_3p2,l_ep2,l_3f2,l_1f3,l_3d3,mold=l_1s0np)

    l_1s0pp = pw_param(:,1)
    l_1s0np = l_1s0pp + pw_param(:,2)
    l_3p0 = pw_param(:, 5)
    l_1p1 = pw_param(:, 6)
    l_3p1 = pw_param(:, 7)
    l_3s1 = pw_param(:, 8)
    l_ep1 = pw_param(:, 9)
    l_3d1 = pw_param(:,10)
    l_1d2 = pw_param(:,11)
    l_3d2 = pw_param(:,12)
    l_3p2 = pw_param(:,13)
    l_ep2 = pw_param(:,14)
    l_3f2 = pw_param(:,15)
    l_1f3 = pw_param(:,16)
    l_3d3 = pw_param(:,18)

    oper_param(:,1) = (1/16._dp)*l_1s0np + (1/8._dp)*l_1s0pp           &
    & + (3/8._dp)*l_3p0 + (3/40._dp)*l_1p1 + (3/16._dp)*l_3s1          &
    & + (3/10._dp)*l_3p2 + (sqrt(6._dp)/5._dp)*l_ep2                   &
    & + (-9/80._dp)*l_3f2 + (-1/80._dp)*l_1f3                           !O_c

    oper_param(:,2) = (1/48._dp)*l_1s0np + (1/24._dp)*l_1s0pp          &
    & + (1/8._dp)*l_3p0 + (-3/40._dp)*l_1p1 + (-3/16._dp)*l_3s1        &
    & + (1/10._dp)*l_3p2 + (sqrt(6._dp)/15._dp)*l_ep2                  &
    & + (-3/80._dp)*l_3f2 + (1/80._dp)*l_1f3                           !O_tau

    oper_param(:,3) = (-1/16._dp)*l_1s0np + (-1/8._dp)*l_1s0pp         &
    & + (1/8._dp)*l_3p0 + (-3/40._dp)*l_1p1 + (1/16._dp)*l_3s1         &
    & + (1/10._dp)*l_3p2 + (sqrt(6._dp)/15._dp)*l_ep2                  &
    & + (-3/80._dp)*l_3f2 + (1/80._dp)*l_1f3                           !O_sigma

    oper_param(:,4) = (-1/48._dp)*l_1s0np + (-1/24._dp)*l_1s0pp        &
    & + (1/24._dp)*l_3p0 + (3/40._dp)*l_1p1 + (-1/16._dp)*l_3s1        &
    & + (1/30._dp)*l_3p2 + (sqrt(6._dp)/45._dp)*l_ep2                  &
    & + (-1/80._dp)*l_3f2 + (-1/80._dp)*l_1f3                           !O_sigma_tau

    oper_param(:,5) = (sqrt(2._dp)/16._dp)*l_ep1                       &
    & + (5*sqrt(6._dp)/48._dp)*l_ep2                                    !O_t

    oper_param(:,6) = (-sqrt(2._dp)/16._dp)*l_ep1                      &
    & + (5*sqrt(6._dp)/144._dp)*l_ep2                                   !O_t_tau

    l_1s0pp = pw_param(:,1)
    l_1s0np = l_1s0pp + pw_param(:,2)
    l_3p0 = pw_param(:, 5)
    l_1p1 = pw_param(:, 6)
    l_3p1 = pw_param(:, 7)
    l_3s1 = pw_param(:, 8)
    l_ep1 = pw_param(:, 9)
    l_3d1 = pw_param(:,10)
    l_1d2 = pw_param(:,11)
    l_3d2 = pw_param(:,12)
    l_3p2 = pw_param(:,13)
    l_ep2 = pw_param(:,14)
    l_3f2 = pw_param(:,15)
    l_1f3 = pw_param(:,16)
    l_3d3 = pw_param(:,18)

    oper_param(:,7) = (-3/8._dp)*l_3p1 + (sqrt(2._dp)/56._dp)*l_ep1    &
    & + (-1/40._dp)*l_3d1 + (-1/24._dp)*l_3d2 + (3/8._dp)*l_3p2        &
    & + (sqrt(6._dp)/8._dp)*l_ep2 + (1/15._dp)*l_3d3                    !O_ls

    oper_param(:,8) = (-1/8._dp)*l_3p1 + (-sqrt(2._dp)/56._dp)*l_ep1   &
    & + (1/40._dp)*l_3d1 + (1/24._dp)*l_3d2  + (1/8._dp)*l_3p2         &
    & + (sqrt(6._dp)/24._dp)*l_ep2 + (-1/15._dp)*l_3d3                  !O_ls_tau

    oper_param(:,9) = (-1/96._dp)*l_1s0np + (-1/48._dp)*l_1s0pp        &
    & + (-9/32._dp)*l_3p0 + (-1/160._dp)*l_1p1 + (9/32._dp)*l_3p1      &
    & + (-1/32._dp)*l_3s1 + (-sqrt(2._dp)/56._dp)*l_ep1                &
    & + (-1/160._dp)*l_3d1 + (1/32._dp)*l_1d2 + (1/32._dp)*l_3d2       &
    & + (-9/160._dp)*l_3p2 + (-9*sqrt(6._dp)/40._dp)*l_ep2             &
    & + (9/160._dp)*l_3f2 + (1/160._dp)*l_1f3 + (1/160._dp)*l_3d3       !O_l2

    oper_param(:,10) = (1/96._dp)*l_1s0np + (1/48._dp)*l_1s0pp         &
    & + (-3/32._dp)*l_3p0 + (1/160._dp)*l_1p1 + (3/32._dp)*l_3p1       &
    & + (-1/96._dp)*l_3s1 + (-sqrt(2._dp)/168._dp)*l_ep1               &
    & + (-1/480._dp)*l_3d1 + (-1/32._dp)*l_1d2 + (1/96._dp)*l_3d2      &
    & + (-3/160._dp)*l_3p2 + (-3*sqrt(6._dp)/40._dp)*l_ep2             &
    & + (3/160._dp)*l_3f2 + (-1/160._dp)*l_1f3 + (1/480._dp)*l_3d3      !O_l2_sigma

    oper_param(:,11) = (-1/288._dp)*l_1s0np + (-1/144._dp)*l_1s0pp     &
    & + (-3/32._dp)*l_3p0 + (1/160._dp)*l_1p1 + (3/32._dp)*l_3p1       &
    & + (1/32._dp)*l_3s1 + (sqrt(2._dp)/56._dp)*l_ep1                  &
    & + (1/160._dp)*l_3d1 + (1/96._dp)*l_1d2 + (-1/32._dp)*l_3d2       &
    & + (-3/160._dp)*l_3p2 + (-3*sqrt(6._dp)/40._dp)*l_ep2             &
    & + (3/160._dp)*l_3f2 + (-1/160._dp)*l_1f3 + (-1/160._dp)*l_3d3     !O_l2_tau

    oper_param(:,12) = (1/288._dp)*l_1s0np + (1/144._dp)*l_1s0pp       &
    & + (-1/32._dp)*l_3p0 + (-1/160._dp)*l_1p1 + (1/32._dp)*l_3p1      &
    & + (1/96._dp)*l_3s1 + (sqrt(2._dp)/168._dp)*l_ep1                 &
    & + (1/480._dp)*l_3d1 + (-1/96._dp)*l_1d2 + (-1/96._dp)*l_3d2      &
    & + (-1/160._dp)*l_3p2 + (-sqrt(6._dp)/40._dp)*l_ep2               &
    & + (1/160._dp)*l_3f2 + (1/160._dp)*l_1f3 + (-1/480._dp)*l_3d3      !O_l2_sigma_tau

    oper_param(:,13) = (1/4._dp)*l_3p0 + (-3/8._dp)*l_3p1              &
    & + (sqrt(2._dp)/28._dp)*l_ep1 + (1/40._dp)*l_3d1                  &
    & + (-1/24._dp)*l_3d2 + (1/8._dp)*l_3p2                            &
    & + (sqrt(6._dp)/4._dp)*l_ep2 + (1/60._dp)*l_3d3                    !O_ls2

    oper_param(:,14) = (1/12._dp)*l_3p0 + (-1/8._dp)*l_3p1             &
    & + (-sqrt(2._dp)/28._dp)*l_ep1 + (-1/40._dp)*l_3d1                &
    & + (1/24._dp)*l_3d2 + (1/24._dp)*l_3p2                            &
    & + (sqrt(6._dp)/12._dp)*l_ep2 + (-1/60._dp)*l_3d3                  !O_ls2_tau

    oper_param(:,15) = (-1/24._dp)*l_1s0np + (1/24._dp)*l_1s0pp         !O_T

    oper_param(:,16) = -oper_param(:,15)                                !O_sigma_T
    oper_param(:,17:19) = 0._dp                                         !O_t_T, O_tau_z, O_sigma_tau_z
    oper_param(:,20) = -oper_param(:,15)/6._dp                          !O_l2_T
    oper_param(:,21) = -oper_param(:,20)                                !O_l2_sigma_T

    oper_param = oper_param/proton_mass_MeV*(hbar_c_MeV_fm**2)          !Returned in units of MeV*fm
end subroutine partial_waves_2_operators

subroutine partial_wave_decomp(oper_param,s,t,j,l,lp,tz1,tz2,lambdas)
    implicit none
    real(dp), intent(in) :: oper_param(:,:)
    integer, intent(in) :: s,t,j,l,lp,tz1,tz2
    real(dp), intent(out) :: lambdas(:)

    integer :: s1ds2,t1dt2,l2,ls,T12, n_lambdas, oper_param_shape(1:2)

    real(dp), allocatable, dimension(:) :: vc,vtau,vsigma,vsigmatau,   &
    & S12,vt,vttau,vls,vlstau,vl2,vl2tau,vl2sigma,vl2sigmatau,vls2,    &
    & vls2tau,viT,viTsigma,vl2iT,vl2iTsigma

    oper_param_shape = shape(oper_param)
    n_lambdas = size(lambdas)

    if (n_lambdas .ne. oper_param_shape(1)) then
        print*, 'Number of lambdas are not the same in partial_wave_decomp'
        stop
    endif

    allocate(vc(1:n_lambdas))
    allocate(vtau,vsigma,vsigmatau,S12,vt,vttau,vls,vlstau,vl2,vl2tau, &
    &        vl2sigma,vl2sigmatau,vls2,vls2tau,viT,viTsigma,vl2iT,     &
    &        vl2iTsigma,mold=vc)

    s1ds2 = 4*s-3
    t1dt2 = 4*t-3
    l2 = delta_k(l,lp)*l*(l+1)
    ls = delta_k(l,lp)*(j*(j+1)-l*(l+1)-s*(s+1))/2
    T12 = delta_k(l,lp)*delta_k(t,1)*(3*tz1*tz2-t1dt2)
    
    S12 = delta_k(s,1)*(2*delta_k(j,l)*delta_k(j,lp)+                  &
        & (-2*(j-1)*delta_k(j,l+1)*delta_k(j,lp+1)                     &
        &  -2*(j+2)*delta_k(j+1,l)*delta_k(j+1,lp)                     &
        &  +6*sqrt(j*(j+1._dp))*( delta_k(j+1,l)*delta_k(j-1,lp)       &
        &                        +delta_k(j-1,l)*delta_k(j+1,lp)))     &
        & /(2*j+1._dp))
    
    vc = oper_param(:,1)
    vtau = oper_param(:,2)
    vsigma = oper_param(:,3)
    vsigmatau = oper_param(:,4)
    vt = oper_param(:,5)
    vttau = oper_param(:,6)
    vls = oper_param(:,7)
    vlstau = oper_param(:,8)
    vl2 = oper_param(:,9)
    vl2sigma = oper_param(:,10)
    vl2tau = oper_param(:,11)
    vl2sigmatau = oper_param(:,12)
    vls2 = oper_param(:,13)
    vls2tau = oper_param(:,14)
    viT = oper_param(:,15)
    viTsigma = oper_param(:,16)
    vl2iT = oper_param(:,20)
    vl2iTsigma = oper_param(:,21)
    
    lambdas = (1-(-1)**(l+s+t))/2._dp *                                 &
    & (delta_k(l,lp)*(vc+t1dt2*vtau+s1ds2*vsigma+t1dt2*s1ds2*vsigmatau) &
    &  +S12*(vt+t1dt2*vttau)+ls*(vls+t1dt2*vlstau)                      &
    &  +l2*(vl2+t1dt2*vl2tau+s1ds2*vl2sigma+t1dt2*s1ds2*vl2sigmatau)    &
    &  +ls**2*(vls2+t1dt2*vls2tau)+T12*(viT+s1ds2*viTsigma)             &
    &  +T12*l2*(vl2iT+s1ds2*vl2iTsigma))

end subroutine partial_wave_decomp

function delta_k(i,j) result(delta)
    implicit none
    integer, intent(in) :: i,j
    integer :: delta

    if (i.eq.j) then
        delta = 1
    else
        delta = 0
    endif

end function delta_k

end module c2m_basis