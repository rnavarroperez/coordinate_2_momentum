program coordinate_2_momentum    

use types
use c2m_parameters
use c2m_io

implicit none

real(dp) :: parameters(1:n_lambdas,1:n_waves)
    
call read_parameters(parameters)
! call partial_waves_2_operators
! call write_momentum_dependence
    
end program coordinate_2_momentum