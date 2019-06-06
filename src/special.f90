module special

! Special functions:
! This module offers some special functions such as
!- spherical Bessel functions
use types
implicit none

private
public spherical_bessel_jn

contains

real(dp) function spherical_bessel_jn(n, x) result(r)
integer, intent(in) :: n
real(dp), intent(in) :: x
integer :: nm
real(dp) :: sj(0:n), dj(0:n)
call sphj(n, x, nm, sj, dj)
if (nm /= n) then
    print*, "spherical_bessel_jn: sphj didn't converge"
    stop
endif
r = sj(n)
end function

! The SPHJ, SPHY, MSTA1, MSTA2 routines below are taken from SciPy's specfun.f.
! Authors: Shanjie Zhang and Jianming Jin
! Copyrighted but permission granted to use code in programs.
! They have been refactored into modern Fortran.
        subroutine sphj(n,x,nm,sj,dj)
!       =======================================================
!       Purpose: Compute spherical Bessel functions jn(x) and
!                their derivatives
!       Input :  x --- Argument of jn(x)
!                n --- Order of jn(x)  ( n = 0,1,… )
!       Output:  SJ(n) --- jn(x)
!                DJ(n) --- jn'(x)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       =======================================================

        integer,intent(in) :: n
        real(dp),intent(in) :: x
        integer,intent(out) :: nm
        real(dp),dimension(0:n),intent(out) :: sj
        real(dp),dimension(0:n),intent(out) :: dj

        integer :: k,m
        real(dp) :: sa,sb,f,f0,f1,cs

        ! .3333333333333333d0 in original
        real(dp),parameter :: one_third = 1.0_dp / 3.0_dp

        nm=n
        if (abs(x)<1.0e-100_dp) then
           do k=0,n
              sj(k)=0.0_dp
              dj(k)=0.0_dp
           end do
           sj(0)=1.0_dp
           if (n>0) then
              dj(1)=one_third
           endif
           return
        endif
        sj(0)=sin(x)/x
        dj(0)=(cos(x)-sin(x)/x)/x
        if (n<1) then
           return
        endif
        sj(1)=(sj(0)-cos(x))/x
        if (n>=2) then
           sa=sj(0)
           sb=sj(1)
           m=msta1(x,200)
           if (m<n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f=0.0_dp
           f0=0.0_dp
           f1=1.0_dp-100
           do k=m,0,-1
              f=(2.0_dp*k+3.0_dp)*f1/x-f0
              if (k<=nm) sj(k)=f
              f0=f1
              f1=f
           end do
           cs=0.0_dp
           if (abs(sa)>abs(sb)) cs=sa/f
           if (abs(sa)<=abs(sb)) cs=sb/f0
           do k=0,nm
              sj(k)=cs*sj(k)
           end do
        endif
        do k=1,nm
           dj(k)=sj(k-1)-(k+1.0_dp)*sj(k)/x
        end do
        end subroutine sphj

        subroutine sphy(n,x,nm,sy,dy)
!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ≥ 0 )
!                n --- Order of yn(x) ( n = 0,1,… )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================

        integer,intent(in) :: n
        real(dp),intent(in) :: x
        integer,intent(out) :: nm
        real(dp),dimension(0:n),intent(out) :: sy
        real(dp),dimension(0:n),intent(out) :: dy

        integer :: k
        real(dp) :: f0,f1,f

        nm=n
        if (x<1.0e-60_dp) then
           do k=0,n
              sy(k)=-1.0e+300_dp
              dy(k)=1.0e+300_dp
           end do
           return
        endif
        sy(0)=-cos(x)/x
        f0=sy(0)
        dy(0)=(sin(x)+cos(x)/x)/x
        if (n<1) then
           return
        endif
        sy(1)=(sy(0)-sin(x))/x
        f1=sy(1)
        do k=2,n
           f=(2.0_dp*k-1.0_dp)*f1/x-f0
           sy(k)=f
           if (abs(f)>=1.0e+300_dp) exit
           f0=f1
           f1=f
        end do
        nm=k-1
        do k=1,nm
           dy(k)=sy(k-1)-(k+1.0_dp)*sy(k)/x
        end do
        end subroutine sphy

        integer function msta1(x,mp)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

        real(dp),intent(in) :: x
        integer,intent(in) :: mp

        real(dp) :: a0,f0,f1,f
        integer :: n0,n1,it,nn

        a0=abs(x)
        n0=int(1.1_dp*a0)+1
        f0=envj(n0,a0)-mp
        n1=n0+5
        f1=envj(n1,a0)-mp
        do it=1,20
           nn=int(n1-(n1-n0)/(1.0_dp-f0/f1))
           f=envj(nn,a0)-mp
           if(abs(nn-n1)<1) exit
           n0=n1
           f0=f1
           n1=nn
           f1=f
        end do
        msta1=nn
        end function msta1

        integer function msta2(x,n,mp)
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

        real(dp),intent(in) :: x
        integer,intent(in) :: n
        integer,intent(in) :: mp

        real(dp) :: a0,hmp,ejn,obj,f0,f1,f
        integer :: n0,n1,nn,it

        a0=abs(x)
        hmp=0.5_dp*mp
        ejn=envj(n,a0)
        if (ejn<=hmp) then
           obj=mp
           n0=int(1.1_dp*a0)+1
        else
           obj=hmp+ejn
           n0=n
        endif
        f0=envj(n0,a0)-obj
        n1=n0+5
        f1=envj(n1,a0)-obj
        do it=1,20
           nn=int(n1-(n1-n0)/(1.0_dp-f0/f1))
           f=envj(nn,a0)-obj
           if (abs(nn-n1)<1) exit
           n0=n1
           f0=f1
           n1=nn
           f1=f
        end do
        msta2=nn+10
        end function msta2

        real(dp) function envj(n, x) result(r)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        r = log10(6.28_dp*n)/2 - n*log10(1.36_dp*x/n)
        end function

end module

