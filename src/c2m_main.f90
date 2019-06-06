program coordinate_2_momentum    

use types
use c2m_io
use c2m_basis

implicit none

integer, parameter :: n_lambdas   = 3 !< number of strength parameters in the parameters array
integer, parameter :: n_waves     = 20!< number of paratial waves in the parameters array
integer, parameter :: n_operators = 21!< number of operators in the potential basis array

real(dp) :: pw_parameters(1:n_lambdas,1:n_waves)
real(dp) :: oper_parameters(1:n_lambdas,1:n_operators)
real(dp), parameter :: delta_r = 0.6_dp
    
call read_parameters(pw_parameters)
call partial_waves_2_operators(pw_parameters,oper_parameters)
call write_momentum_dependence(oper_parameters,delta_r)
    
end program coordinate_2_momentum