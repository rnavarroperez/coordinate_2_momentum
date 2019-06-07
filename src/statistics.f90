module statistics
use types, only: dp
implicit none

private
public mean, sample_variance

contains

real(dp) function mean(x) result(mu)
    implicit none
    real(dp), intent(in) :: x(:)
    if (size(x).lt.1) then
        print*, 'Array size needs to be at least 1 in mean'
        stop
    endif
    mu = sum(x)/size(x)
end function mean

real(dp) function sample_variance(x) result(var)
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp) :: mu
    if (size(x).lt.2) then
        print*, 'Array size needs to be at least 2 in sample_variance'
        stop
    endif
    mu = mean(x)
    var = sum((mu-x)**2)/(size(x)-1)
end function sample_variance
    
end module statistics