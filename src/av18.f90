module av18
use types, only: dp 
use constants, only: hc=>hbar_c_MeV_fm, mpic=>pion_c_mass_MeV, mpi0=>pion_0_mass_MeV, mpi=>pion_mass_MeV
implicit none
private
public av18_oper_basis, n_operators

integer, parameter :: n_parameters = 44
real(dp), parameter :: f2=0.075_dp !< pion nucleon coupling constant
real(dp), parameter, dimension(1:n_parameters) :: av18_parameters=[&
    &     -7.62701, &
    &   1815.49200, &
    &   1847.80590, &
    &   1813.53150, &
    &   1811.57100, &
    &      1.07985, &
    &   -190.09490, &
    &   -811.20400, &
    &     -0.62697, &
    &   -570.55710, &
    &    819.12220, &
    &      0.06709, &
    &    342.06690, &
    &   -615.23390, &
    &      0.74129, &
    &      9.34180, &
    &   -376.43840, &
    &     -8.62770, &
    &   2605.26820, &
    &    441.97330, &
    &      1.48560, &
    &  -1126.83590, &
    &    370.13240, &
    &      0.10180, &
    &     86.06580, &
    &   -356.51750, &
    &     -0.13201, &
    &    253.43500, &
    &     -1.00760, &
    &      0.07357, &
    &   -217.57910, &
    &     18.39350, &
    &    -11.27028, &
    &   3346.68740, &
    &    -10.66788, &
    &   3126.55420, &
    &    -11.27028, &
    &   3342.76640, &
    &      0.12472, &
    &     16.77800, &
    &     -2.09971, &
    &   1204.43010, &
    &     -0.31452, &
    &    217.45590   ] 

integer, parameter :: n_operators=18

contains

subroutine av18_oper_basis(r,v_operator)
    implicit none
    real(dp), intent(in) ::  r
    real(dp), intent(out) :: v_operator(:)

    real(dp), parameter :: &
        mu0 = mpi0/hc, &
        muc = mpic/hc, &
        mu  = mpi/hc,  &
        cpi = 2.1_dp, &
        rws = 0.5_dp, &
        aiws = 5._dp, &
        small = 1.0e-4_dp

    real(dp) :: x,x0,xc,ypi,tpi,ypi0,tpi0,ypic,tpic,rcut,ws,ws0,wsp,&
        & wsx,wsx2,dypi00,dypic0,ypi0p,ypicp,tpi2,p11pp,p11np,p11nn,&
        & pt1pp,pt1np,pt1nn,pls1,pl211,pls21,p10,pt0,pls0,pl210,pls20,&
        & p01pp,p01np,p01nn,pl201,p00,pl200,p11,p11cd,p11cs,pt1,pt1cd,&
        & pt1cs,p01,p01cd,p01cs

    if(size(v_operator).ne.n_operators) then
        print*, 'v_operator needs to be size 18 in av18_oper_basis'
        stop
    endif

    v_operator = 0._dp
    x = mu*r
    x0 = mu0*r
    xc = muc*r
    if (r.le.small) then
        tpi=3*cpi**2*r/mu**3
        ypi0=(mpi0/mpic)**2*(mpi0/3)*cpi*r/mu0
        tpi0=3*cpi*ypi0/mu0**2
        ypic=(mpic/3)*cpi*r/muc
        tpic=3*cpi*ypic/muc**2
    else
        rcut=1-exp(-cpi*r*r) 
        ypi=exp(-x)*rcut/x   
        tpi=(1+(3+3/x)/x)*ypi*rcut 
        ypi0=(mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        tpi0=(1+(3+3/x0)/x0)*ypi0*rcut
        ypic=(mpic/3)*exp(-xc)*rcut/xc
        tpic=(1+(3+3/xc)/xc)*ypic*rcut
    end if
    ypi0=f2*ypi0
    ypic=f2*ypic
    tpi0=f2*tpi0
    tpic=f2*tpic
    tpi2=tpi*tpi 
    ws=1/(1+exp((r-rws)*aiws)) 
    ws0=1/(1+exp(-rws*aiws))   
    wsp=ws*(1+aiws*exp(-rws*aiws)*ws0*r) 
    wsx=ws*x
    wsx2=wsx*x
    dypi00=(mpi0/mpic)**2*(mpi0/3)*cpi/mu0 
    dypic0=(mpic/3)*cpi/muc
    ypi0p=ypi0-f2*dypi00*ws*r/ws0
    ypicp=ypic-f2*dypic0*ws*r/ws0
    ypi=(ypi0+2*ypic)/3
    tpi=(tpi0+2*tpic)/3

    p11pp=  av18_parameters( 1)*tpi2 + av18_parameters( 2)*wsp + av18_parameters( 3)*wsx2 + ypi0p
    p11np=  av18_parameters( 1)*tpi2 + av18_parameters( 4)*wsp + av18_parameters( 3)*wsx2 - ypi0p+2*ypicp
    p11nn=  av18_parameters( 1)*tpi2 + av18_parameters( 5)*wsp + av18_parameters( 3)*wsx2 + ypi0p
    pt1pp=  av18_parameters( 6)*tpi2 + av18_parameters( 7)*wsx + av18_parameters( 8)*wsx2 + tpi0
    pt1np=  av18_parameters( 6)*tpi2 + av18_parameters( 7)*wsx + av18_parameters( 8)*wsx2 - tpi0+2*tpic
    pt1nn=  av18_parameters( 6)*tpi2 + av18_parameters( 7)*wsx + av18_parameters( 8)*wsx2 + tpi0
    pls1=   av18_parameters( 9)*tpi2 + av18_parameters(10)*wsp + av18_parameters(11)*wsx2
    pl211=  av18_parameters(12)*tpi2 + av18_parameters(13)*wsp + av18_parameters(14)*wsx2
    pls21=  av18_parameters(15)*tpi2 + av18_parameters(16)*wsp + av18_parameters(17)*wsx2
    p10=    av18_parameters(18)*tpi2 + av18_parameters(19)*wsp + av18_parameters(20)*wsx2 - ypi0p-2*ypicp
    pt0=    av18_parameters(21)*tpi2 + av18_parameters(22)*wsx + av18_parameters(23)*wsx2 - tpi0-2*tpic
    pls0=   av18_parameters(24)*tpi2 + av18_parameters(25)*wsp + av18_parameters(26)*wsx2
    pl210=  av18_parameters(27)*tpi2 + av18_parameters(28)*wsp + av18_parameters(29)*wsx2
    pls20=  av18_parameters(30)*tpi2 + av18_parameters(31)*wsp + av18_parameters(32)*wsx2
    p01pp=  av18_parameters(33)*tpi2 + av18_parameters(34)*wsp - 3*ypi0p
    p01np=  av18_parameters(35)*tpi2 + av18_parameters(36)*wsp - 3*(-ypi0p+2*ypicp)
    p01nn=  av18_parameters(37)*tpi2 + av18_parameters(38)*wsp - 3*ypi0p
    pl201=  av18_parameters(39)*tpi2 + av18_parameters(40)*wsp
    p00=    av18_parameters(41)*tpi2 + av18_parameters(42)*wsp - 3*(-ypi0p-2*ypicp)
    pl200=  av18_parameters(43)*tpi2 + av18_parameters(44)*wsp

    ! p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p
    ! p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp
    ! p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p
    ! pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
    ! pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic
    ! pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
    ! pls1=    -.62697*tpi2 -570.5571*wsp +819.1222*wsx2
    ! pl211=    .06709*tpi2 +342.0669*wsp -615.2339*wsx2
    ! pls21=    .74129*tpi2   +9.3418*wsp -376.4384*wsx2
    ! p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp
    ! pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic
    ! pls0=     .10180*tpi2  +86.0658*wsp -356.5175*wsx2
    ! pl210=   -.13201*tpi2 +253.4350*wsp   -1.0076*wsx2
    ! pls20=    .07357*tpi2 -217.5791*wsp  +18.3935*wsx2
    ! p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p
    ! p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
    ! p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p
    ! pl201=    .12472*tpi2  +16.7780*wsp
    ! p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp)
    ! pl200=   -.31452*tpi2 +217.4559*wsp

    p11=(p11pp+p11nn+p11np)/3
    p11cd=(.5_dp*(p11pp+p11nn)-p11np)/6
    p11cs=(p11pp-p11nn)/4
    pt1=(pt1pp+pt1nn+pt1np)/3
    pt1cd=(.5_dp*(pt1pp+pt1nn)-pt1np)/6
    pt1cs=(pt1pp-pt1nn)/4
    p01=(p01pp+p01nn+p01np)/3
    p01cd=(.5_dp*(p01pp+p01nn)-p01np)/6
    p01cs=(p01pp-p01nn)/4

    v_operator(1)=.0625_dp*(9*p11+3*p10+3*p01+p00)
    v_operator(2)=.0625_dp*(3*p11-3*p10  +p01-p00)
    v_operator(3)=.0625_dp*(3*p11  +p10-3*p01-p00)
    v_operator(4)=.0625_dp*(  p11  -p10  -p01+p00)
    v_operator(5)=.25_dp*(3*pt1+pt0)
    v_operator(6)=.25_dp*(  pt1-pt0)
    v_operator(7)=.25_dp*(3*pls1+pls0)
    v_operator(8)=.25_dp*(  pls1-pls0)
    v_operator(9)= .0625_dp*(9*pl211+3*pl210+3*pl201+pl200)
    v_operator(10)=.0625_dp*(3*pl211-3*pl210+  pl201-pl200)
    v_operator(11)=.0625_dp*(3*pl211+  pl210-3*pl201-pl200)
    v_operator(12)=.0625_dp*(  pl211-  pl210-  pl201+pl200)
    v_operator(13)=.25_dp*(3*pls21+pls20)
    v_operator(14)=.25_dp*(  pls21-pls20)
    v_operator(15)=.25_dp*(3*p11cd+p01cd)
    v_operator(16)=.25_dp*(  p11cd-p01cd)
    v_operator(17)=pt1cd
    v_operator(18)=p01cs
end subroutine av18_oper_basis

subroutine av18_spin_isospin_basis(tz1,tz2,s,t,v_oper,v_c,v_t,v_ls,v_l2,v_ls2)
    implicit none
    integer, intent(in) :: tz1,tz2,s,t
    real(dp), intent(in) :: v_oper(:)
    real(dp), intent(out) :: v_c,v_t,v_ls,v_l2,v_ls2
    integer :: s1ds2, t1dt2,t12
    s1ds2=4*s-3
    t1dt2=4*t-3
    t12=3*tz1*tz2-t1dt2
    v_c=v_oper(1)+t1dt2*v_oper(2)+s1ds2*v_oper(3)+s1ds2*t1dt2*v_oper(4)&
        &+t12*v_oper(15)+s1ds2*t12*v_oper(16)+(tz1+tz2)*v_oper(18)
    v_t=v_oper(5)+t1dt2*v_oper(6)+t12*v_oper(17)
    v_ls=v_oper(7)+t1dt2*v_oper(8)
    v_l2=v_oper(9)+t1dt2*v_oper(10)+s1ds2*v_oper(11)+s1ds2*t1dt2*v_oper(12)
    v_ls2=v_oper(13)+t1dt2*v_oper(14)
end subroutine av18_spin_isospin_basis

end module av18