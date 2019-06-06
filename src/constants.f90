module constants
use types, only: dp
implicit none
private
public pi, e_, i_,hbar_c_MeV_fm,proton_mass_MeV,neutron_mass_MeV

! Constants contain more digits than double precision, so that
! they are rounded correctly. Single letter constants contain underscore so
! that they do not clash with user variables ("e" and "i" are frequently used as
! loop variables)
real(dp), parameter :: pi    = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_    = 2.7182818284590452353602874713527_dp
complex(dp), parameter :: i_ = (0, 1)

! Values from National Institute of Standards and Technology (NIST)
! https://www.nist.gov/pml/fundamental-physical-constants
! Retrieved on June 6th 2019
real(dp), parameter :: hbar_c_MeV_fm    = 197.3269804_dp  !< hbar c in units of MeVfm
real(dp), parameter :: proton_mass_MeV  = 938.27208816_dp !< proton mass in units of MeV
real(dp), parameter :: neutron_mass_MeV = 939.56542052_dp   !< neutron mass in units of MeV
end module
