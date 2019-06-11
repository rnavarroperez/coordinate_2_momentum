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