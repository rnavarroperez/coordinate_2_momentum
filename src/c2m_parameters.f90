module c2m_parameters
use types
implicit none
integer, parameter :: n_waves   = 20!< number of paratial waves in the parameters array
integer, parameter :: n_lambdas = 3!< number of lambda strength parameters in the parameters array

character(len=25), parameter :: f_name = 'tpe218pwanew.dat' !< file nane subfix
  
end module c2m_parameters