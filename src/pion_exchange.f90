module pion_exchange
use types
use constants, only: hbar_c_MeV_fm, pion_c_mass_MeV, pion_0_mass_MeV, pion_mass_MeV,pi,proton_mass_MeV,neutron_mass_MeV
use av18, only: n_operators
use special, only: bessel_k0,bessel_k1
implicit none
private
public v_one_pion_exch, v_two_pion_exch_nlo, v_two_pion_exch_n2lo
real(dp), parameter :: f2pp =  0.075_dp,f2nn =  0.075_dp, f2np = -0.075_dp, f2c  =  0.075_dp
real(dp), parameter :: hc=hbar_c_MeV_fm, mpic = pion_c_mass_MeV, mpi0=pion_0_mass_MeV, mpi=pion_mass_MeV
real(dp), parameter :: mp=proton_mass_MeV,mn=neutron_mass_MeV
real(dp), parameter :: F_pi = 185._dp !< Pion decay constant, value from PhysRevLett.82.4992
contains

subroutine v_two_pion_exch_n2lo(r,lecs,v_oper_basis)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(in) :: lecs(1:3)
    real(dp), intent(out) :: v_oper_basis(1:n_operators)

    real(dp) :: ap_v(1:6,1:4), ap_w(1:6,1:4)
    real(dp) :: v_n2lo(1:4),w_n2lo(1:4)
    real(dp) :: f,x,Xi

    v_oper_basis = 0._dp

    x = mpi*r/hc
    Xi = mpi/mpic

    call construct_ap_tables(lecs,ap_v,ap_w)
    call sum_n2lo_pot(x,ap_v,v_n2lo)
    call sum_n2lo_pot(x,ap_w,w_n2lo)

    f = sqrt(f2c)
    v_n2lo  =  v_n2lo*f**4*Xi**4*mpi
    w_n2lo  =  w_n2lo*f**4*Xi**4*mpi

    v_oper_basis(1) = v_n2lo(1)
    v_oper_basis(2) = w_n2lo(1)
    v_oper_basis(3) = v_n2lo(2)
    v_oper_basis(4) = w_n2lo(2)
    v_oper_basis(5) = v_n2lo(3)
    v_oper_basis(6) = w_n2lo(3)
    v_oper_basis(7) = v_n2lo(4)
    v_oper_basis(8) = w_n2lo(4)
end subroutine v_two_pion_exch_n2lo

subroutine construct_ap_tables(lecs,apv,apw)
    implicit none
    real(dp), intent(in) ::  lecs(1:3)
    real(dp), intent(out) :: apv(1:6,1:4),apw(1:6,1:4)
    real(dp) :: f,gAt,c0t,c1t,c3t,c4t,c04t,m,dgAt

    apv(:,1) = [0.75_dp, 9.0_dp, 27.00_dp, 49.5_dp, 54._dp, 27._dp]
    apv(:,2) = [0.00_dp,-3.0_dp, -9.00_dp,-16.5_dp,-18._dp, -9._dp]
    apv(:,3) = [0.00_dp, 1.5_dp,  6.75_dp, 15.0_dp, 18._dp,  9._dp]
    apv(:,4) = [0.00_dp, 0.0_dp,-12.00_dp,-36.0_dp,-48._dp,-24._dp]

    apw(:,1) = [4.50_dp,12.0_dp, 42.00_dp, 93.0_dp,108._dp, 54._dp]
    apw(:,2) = [0.00_dp,-2.0_dp,-14.00_dp,-31.0_dp,-36._dp,-18._dp]
    apw(:,3) = [0.00_dp, 1.0_dp,  8.50_dp, 26.0_dp, 36._dp, 18._dp]
    apw(:,4) = [0.00_dp, 0.0_dp,  0.00_dp, 24.0_dp, 48._dp, 24._dp]

    m = (mp+mn)/2

    f = sqrt(f2c)
    gAt = F_pi*sqrt(4*pi)*f/mpic

    c0t = 1/gAt**2
    c1t = lecs(1)*m/gAt**2
    c3t = lecs(2)*m/gAt**2
    c4t = lecs(3)*m/gAt**2
    c04t = c0t + 4*c4t

    apw = apw/3.d0

    apv(2,1) = apv(2,1) + 48*c1t +  24*c3t
    apv(3,1) = apv(3,1) + 96*c1t +  96*c3t
    apv(4,1) = apv(4,1) + 48*c1t + 240*c3t
    apv(5,1) = apv(5,1) +          288*c3t
    apv(6,1) = apv(6,1) +          144*c3t

    apw(2,1) = apw(2,1)-  2*c0t
    apw(3,1) = apw(3,1)-  8*c0t
    apw(4,1) = apw(4,1)- 20*c0t
    apw(5,1) = apw(5,1)- 24*c0t
    apw(6,1) = apw(6,1)- 12*c0t
    apw(3,2) = apw(3,2)+  8*c04t/3.d0
    apw(4,2) = apw(4,2)+ 20*c04t/3.d0
    apw(5,2) = apw(5,2)+  8*c04t
    apw(6,2) = apw(6,2)+  4*c04t
    apw(3,3) = apw(3,3)-  4*c04t/3.d0
    apw(4,3) = apw(4,3)- 16*c04t/3.d0
    apw(5,3) = apw(5,3)-  8*c04t
    apw(6,3) = apw(6,3)-  4*c04t
    apw(4,4) = apw(4,4)-  8*c0t
    apw(5,4) = apw(5,4)- 16*c0t
    apw(6,4) = apw(6,4)-  8*c0t

    apv = apv*mpi/m
    apw = apw*mpi/m
end subroutine construct_ap_tables


subroutine sum_n2lo_pot(x,ap,v)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in) :: ap(1:6,1:4)
    real(dp), intent(out) :: v(1:4)

    integer :: i,j
    v = 0._dp
    do i = 1,4
        do j = 1,6
            v(i) = v(i) + ap(j,i)/(x**j)
        enddo
    enddo
    v = v*exp(-2._dp*x)
end subroutine sum_n2lo_pot

subroutine v_two_pion_exch_nlo(r,v_oper_basis)
    implicit none
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: v_oper_basis(1:n_operators)
    real(dp) :: x,Xi,f,gAt,c0t,k0,k1,vs,vt,wc

    v_oper_basis = 0._dp

    x = mpi*r/hc
    Xi = mpi/mpic
    f = sqrt(f2c)
    gAt = F_pi*sqrt(4*pi)*f/mpic
    c0t = 1._dp/gAt**2
    k0 = bessel_k0(2*x)
    k1 = bessel_k1(2*x)
    vs = 12*k0/x**3 + (12+8*x**2)*k1/x**4
    vt =-12*k0/x**3 - (15+4*x**2)*k1/x**4
    wc = (c0t**2+10*c0t-23-4*x**2) *k0/x**3 &
       & +(c0t**2+10*c0t-23+(4*c0t-12)*x**2)*k1/x**4
    vs = 2*vs/pi*f**4*Xi**4*mpi
    vt = 2*vt/pi*f**4*Xi**4*mpi
    wc = 2*wc/pi*f**4*Xi**4*mpi

    v_oper_basis(2) = wc
    v_oper_basis(3) = vs
    v_oper_basis(5) = vt
end subroutine v_two_pion_exch_nlo

subroutine v_one_pion_exch(r,v_oper_basis)
    implicit none
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: v_oper_basis(1:n_operators)
    real(dp) :: vc_11pp, vc_11np, vc_11nn, vc_01pp, vc_01np, vc_01nn
    real(dp) :: vt_11pp, vt_11np, vt_11nn, vt_01pp, vt_01np, vt_01nn
    real(dp) :: vc_10np, vc_00np, vt_10np, vt_00np

    real(dp) :: vc_11_ci, vc_11_cd, vc_11_ca
    real(dp) :: vc_01_ci, vc_01_cd, vc_01_ca
    real(dp) :: vc_10_ci
    real(dp) :: vc_00_ci

    real(dp) :: vt_11_ci, vt_11_cd
    real(dp) :: vt_10_ci

    v_oper_basis = 0._dp

    call v_ope_st(r,1,1, 1, 1,vc_11pp,vt_11pp)
    call v_ope_st(r,1,1,-1, 1,vc_11np,vt_11np)
    call v_ope_st(r,1,1,-1,-1,vc_11nn,vt_11nn)

    call v_ope_st(r,0,1, 1, 1,vc_01pp,vt_01pp)
    call v_ope_st(r,0,1,-1, 1,vc_01np,vt_01np)
    call v_ope_st(r,0,1,-1,-1,vc_01nn,vt_01nn)

    call v_ope_st(r,1,0,-1, 1,vc_10np,vt_10np)
    call v_ope_st(r,0,0,-1, 1,vc_00np,vt_00np)

    ! Transformation betweeh the ST and operator basis can be found 
    ! in section IV of PhysRevC.51.38

    vc_11_ci = (vc_11pp + vc_11np + vc_11nn)/3
    vc_01_ci = (vc_01pp + vc_01np + vc_01nn)/3
    vc_10_ci = vc_10np
    vc_00_ci = vc_00np

    vt_11_ci = (vt_11pp + vt_11np + vt_11nn)/3
    vt_10_ci = vt_10np
    
    vc_11_cd = ((vc_11pp + vc_11nn)/2 - vc_11np)/6
    vc_01_cd = ((vc_01pp + vc_01nn)/2 - vc_01np)/6

    vt_11_cd = ((vt_11pp + vt_11nn)/2 - vt_11np)/6

    vc_11_ca = (vc_11pp - vc_11nn)/4
    vc_01_ca = (vc_01pp - vc_01nn)/4

    v_oper_basis( 1) = (9*vc_11_ci + 3*vc_10_ci + 3*vc_01_ci + vc_00_ci)/16
    v_oper_basis( 2) = (3*vc_11_ci - 3*vc_10_ci +   vc_01_ci - vc_00_ci)/16
    v_oper_basis( 3) = (3*vc_11_ci +   vc_10_ci - 3*vc_01_ci - vc_00_ci)/16
    v_oper_basis( 4) = (  vc_11_ci -   vc_10_ci -   vc_01_ci + vc_00_ci)/16
    v_oper_basis( 5) = (3*vt_11_ci + vt_10_ci)/4
    v_oper_basis( 6) = (  vt_11_ci - vt_10_ci)/4
    v_oper_basis(15) = (3*vc_11_cd + vc_01_cd)/4
    v_oper_basis(16) = (  vc_11_cd - vc_01_cd)/4
    v_oper_basis(17) = vt_11_cd
    v_oper_basis(18) = (3*vc_11_ca + vc_01_ca)/4

end subroutine v_one_pion_exch

subroutine v_ope_st(r,s,t,tz1,tz2,vc,vt)
    implicit none
    real(dp), intent(in) :: r
    integer, intent(in) :: s,t,tz1,tz2
    real(dp), intent(out) :: vc,vt

    integer :: s1ds2
    real(dp) :: ws0,wt0,ws,wt,wsc,wtc

    s1ds2 = 4*s-3
    ws0 = wsope(mpi0,r)
    wt0 = wtope(mpi0,r)
    if(tz1*tz2.eq.1) then
        if(tz1.eq.1) then !pp
            ws = f2pp*ws0 
            wt = f2pp*wt0 
        else !nn
            ws = f2nn*ws0 
            wt = f2nn*wt0
        endif
    elseif(tz1*tz2.eq.-1) then
        !np channel, break symmetry
        wsc = wsope(mpic,r)
        wtc = wtope(mpic,r)
        ws = f2np*ws0 + (-1)**(t+1)*2*f2c*wsc
        wt = f2np*wt0 + (-1)**(t+1)*2*f2c*wtc
    else
        !wrong tz1 and tz2
        print*, 'wrong NN channel in vopepw'
        stop
    endif
    vc  = s1ds2*ws
    vt  = wt
end subroutine v_ope_st

real(dp) function wsope(m,r) result(ws)
    implicit none
    real(dp), intent(in) :: m,r
    real(dp) ::  x,Xi
    x = m*r/hc
    Xi = m/mpic
    ws = exp(-x)/3._dp/x
    ws = ws*Xi**2*m
end function wsope

real(dp) function wtope(m,r) result(wt)
    implicit none
    real(dp), intent(in) :: m,r
    real(dp) ::  x,Xi
    x = m*r/hc
    Xi = m/mpic
    wt = (1._dp+x+x**2/3._dp)*exp(-x)/x**3
    wt = wt*Xi**2*m
end function wtope
    
end module pion_exchange