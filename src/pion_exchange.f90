module pion_exchange
use types
use constants, only: hbar_c_MeV_fm, pion_c_mass_MeV, pion_0_mass_MeV, pion_mass_MeV
implicit none
private
real(dp), parameter :: f2pp =  0.075_dp,f2nn =  0.075_dp, f2np = -0.075_dp, f2c  =  0.075_dp
real(dp), parameter :: hc=hbar_c_MeV_fm, mpic = pion_c_mass_MeV, mpi0=pion_0_mass_MeV, mpi=pion_mass_MeV
contains

subroutine v_one_pion_exch(r,v_oper_basis)
    implicit none
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: v_oper_basis(:)
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