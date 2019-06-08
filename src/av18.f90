module av18
use types, only: dp 
use constants, only: hbar_c_MeV_fm, pion_c_mass_MeV, pion_0_mass_MeV, pion_mass_MeV
implicit none
private

real(dp), parameter :: hc=hbar_c_MeV_fm, mpic = pion_c_mass_MeV, mpi0=pion_0_mass_MeV, mpi=pion_mass_MeV
real(dp), parameter, dimension(1:44) :: av18_parameters=[&
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

end subroutine av18_oper_basis
    
end module av18