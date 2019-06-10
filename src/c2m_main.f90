program coordinate_2_momentum    

use types, only: dp
use c2m_io, only: read_parameters, read_mc_samples, write_momentum_dependence, write_av18_momentum
use c2m_basis, only: set_fit_flags

implicit none

integer, parameter :: n_lambdas   = 3   !< number of strength parameters in the parameters array
integer, parameter :: n_waves     = 20  !< number of paratial waves in the parameters array
integer, parameter :: n_samples   = 999!< number of monte carlo samples 

real(dp) :: pw_parameters(1:n_lambdas,1:n_waves)
real(dp) :: pw_samples(1:n_lambdas,1:n_waves,1:n_samples)
real(dp), parameter :: delta_r_ds = 0.6_dp !in fm
real(dp), parameter :: delta_r_av18 = 1.0e-2_dp !in fm
real(dp), parameter :: r_max = 20.0_dp !in fm

logical :: fit_flags(1:n_lambdas,1:n_waves)

call read_parameters(pw_parameters)
call set_fit_flags(pw_parameters,fit_flags)
call read_mc_samples(pw_parameters,fit_flags,pw_samples)
call write_momentum_dependence(pw_samples,delta_r_ds)

call write_av18_momentum(delta_r_av18,r_max)
    
end program coordinate_2_momentum