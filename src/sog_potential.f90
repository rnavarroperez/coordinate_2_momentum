module sog_potential
use types, only: dp 
use constants, only: hc=>hbar_c_MeV_fm, mpic=>pion_c_mass_MeV, mpi0=>pion_0_mass_MeV, mpi=>pion_mass_MeV
implicit none
private
public sog_potential_oper_basis, n_operators

integer, parameter :: n_operators = 21
integer, parameter :: n_gaussians = 4

real(dp), parameter :: potential_parameters(1:n_operators,1:n_gaussians) = &
    & transpose(reshape([-19.28330126_dp,   126.28715008_dp,  -648.61345585_dp,   694.49367435_dp, & ! c
                        &  2.36233395_dp,   -25.47505195_dp,   130.03224633_dp,  -284.71844492_dp, & ! tau
                        &  6.05581487_dp,   -75.18832503_dp,   372.41961972_dp,  -530.80008401_dp, & ! sigma
                        &  7.36008330_dp,   -48.55160272_dp,   273.71591816_dp,  -349.00547346_dp, & ! sigma tau
                        &  1.99828652_dp,   -22.12164190_dp,    70.84584496_dp,   -50.72248959_dp, & ! t
                        & 15.02271531_dp,   -38.34776035_dp,   183.80564790_dp,  -160.48060286_dp, & ! t tau
                        & -2.61725312_dp,    39.43014573_dp,  -217.03270342_dp,  -109.64162556_dp, & ! ls
                        &  0.01009424_dp,     2.59116238_dp,   -26.57555840_dp,   -77.56809604_dp, & ! ls tau
                        &  1.43519736_dp,   -23.58906341_dp,    67.86552330_dp,   144.11773134_dp, & ! l2
                        & -0.41138176_dp,     8.33346137_dp,   -82.98819447_dp,   175.09618737_dp, & ! l2 tau
                        & -0.09972181_dp,     2.25339230_dp,   -51.87819771_dp,   175.08890636_dp, & ! l2 sigma
                        & -0.26615545_dp,     6.63257735_dp,   -55.34306118_dp,   100.71528331_dp, & ! l2 sigma tau
                        &  0.46072934_dp,   -11.65544792_dp,   150.58275714_dp,  -302.07573779_dp, & ! ls2
                        &  0.71538487_dp,   -18.88652666_dp,   141.73160452_dp,  -182.73368764_dp, & ! ls2 tau
                        &  0.63788724_dp,    -7.90421846_dp,    24.23180376_dp,   -19.73899169_dp, & ! T
                        & -0.63788724_dp,     7.90421846_dp,   -24.23180376_dp,    19.73899169_dp, & ! sigma T
                        &  0.00000000_dp,     0.00000000_dp,     0.00000000_dp,     0.00000000_dp, & ! tT
                        &  0.00000000_dp,     0.00000000_dp,     0.00000000_dp,     0.00000000_dp, & ! tau_z
                        &  0.00000000_dp,     0.00000000_dp,     0.00000000_dp,     0.00000000_dp, & ! sigma tau_z
                        & -0.10631454_dp,     1.31736974_dp,    -4.03863396_dp,     3.28983195_dp, & ! l2 T
                        &  0.10631454_dp,    -1.31736974_dp,     4.03863396_dp,    -3.28983195_dp  & ! l2 sigma T
                        & ], [n_gaussians,n_operators]))

real(dp), parameter :: range = 2.30347728_dp

real(dp), parameter :: f2=0.075_dp !< pion nucleon coupling constant

contains

subroutine sog_potential_oper_basis(r,v_operators)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out) :: v_operators(:)

    real(dp) :: gaussian 
    integer :: i

    if(size(v_operators).ne.n_operators) then
        print*, 'v_operators needs to be size 21 in sog_potential_oper_basis'
        stop
    endif

    if (r .le. 3._dp) then
        v_operators = 0._dp
        do i=1,n_gaussians
            gaussian = exp(-0.5_dp*((i+1)*r/range)**2)
            v_operators = v_operators + potential_parameters(:,i)*gaussian
        enddo
    else
        call ope_tail(r,v_operators)
    endif
    
end subroutine sog_potential_oper_basis

subroutine ope_tail(r, v_operators)
    implicit none
    real(dp), intent(in) ::  r
    real(dp), intent(out) :: v_operators(:)

    real(dp) :: mu0, muc, ypi0, ypic, tpi0, tpic

    if(size(v_operators).ne.n_operators) then
        print*, 'v_operators needs to be size 21 in ope_tail'
        stop
    endif

    v_operators = 0._dp

    mu0 = mpi0/hc
    muc = mpic/hc
    ypi0 = f2*(mpi0/mpic)**2*mpi0/3._dp*exp(-mu0*r)/(mu0*r)
    ypic = f2*mpic/3._dp*exp(-muc*r)/(muc*r)
    tpi0 = ypi0*(1 + 3/(mu0*r) + 3/(mu0*r)**2)
    tpic = ypic*(1 + 3/(muc*r) + 3/(muc*r)**2)
    v_operators( 4) = (ypi0 + 2*ypic)/3._dp
    v_operators( 6) = (tpi0 + 2*tpic)/3._dp
    v_operators(16) = (ypi0 -   ypic)/3._dp
    v_operators(17) = (tpi0 -   tpic)/3._dp

end subroutine ope_tail
    
end module sog_potential