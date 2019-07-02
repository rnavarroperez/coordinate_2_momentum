module c2m_io

use types, only: dp
use c2m_transform, only: delta_shell_2_momentum,sample_av18,transform_all_oper,n_q_operators, n_av18q_operators
use c2m_montecarlo, only: mc_momentum_dependence
use constants, only: pi
implicit none

private
public :: read_parameters, write_momentum_dependence, read_mc_samples, write_av18_momentum

character(len=25), parameter :: f_name = 'tpe218pwanew.dat' !< file nane subfix

contains

subroutine read_mc_samples(parameters,flags,samples)
    implicit none
    real(dp), intent(in) :: parameters(:,:)
    logical, intent(in) :: flags(:,:)
    real(dp), intent(out) :: samples(:,:,:)

    character(len=128) :: filename
    integer :: samples_shape(1:3), n_lambdas, n_waves, n_samples !array sizes
    integer :: unit, ierr ! file manipulation
    logical :: exists !file manipulation
    integer :: n_active 
    real(dp), allocatable :: mc_lambdas(:)
    integer :: id, pp_data, np_data
    real(dp) :: chi2_pp, chi2_np

    integer :: is,il,iw,iac
    

    samples_shape = shape(samples)
    n_lambdas = samples_shape(1)
    n_waves = samples_shape(2)
    n_samples = samples_shape(3)

    if(any(shape(parameters).ne.shape(flags))) then
        print*, 'Parameters and flags have inconsistent shapes in read_mc_samples'
        stop
    else if ( any( samples_shape(1:2) .ne. shape(flags)) ) then
        print*, 'Samples and flags have inconsistent shapes in read_mc_samples'
        stop        
    endif

    samples = 0._dp

    n_active = count(flags)
    allocate(mc_lambdas(1:n_active))
    
    filename = 'mc_lambdas_fixedZ_0001_1000'//f_name
    inquire(file=trim(filename), exist=exists)
    if (exists) then
        open(newunit = unit, file=trim(filename), action="read", iostat=ierr)
        if (ierr .eq. 0) then
            read(unit,*) !skipping first line (this contains the central values in the mc sample)
            do is=1,n_samples
                read(unit,*) id,pp_data,np_data,chi2_pp,chi2_np,mc_lambdas(:)
                iac = 0
                do iw=1,n_waves
                    do il=1,n_lambdas
                        if (flags(il,iw)) then
                            iac = iac + 1
                            samples(il,iw,is) = mc_lambdas(iac)
                        endif
                    enddo
                enddo
            enddo
        else
            print*, "Error ", ierr ," attempting to open file ", trim(filename)
            stop
        endif
        close(unit)
    else
        print*, "Error -- cannot find file: ", trim(filename)
        stop
    endif

end subroutine read_mc_samples

subroutine read_parameters(param)
    implicit none
    real(dp), intent(out) :: param(:,:) !< array to be filled with the potential parameters in partial wave basis

    integer            :: unit, ierr,readerr,i_waves,i_lambdas
    character(len=128) :: filename
    logical            :: exists
    integer            :: param_shape(1:2),n_lambdas,n_waves

    param_shape = shape(param)
    n_lambdas = param_shape(1)
    n_waves = param_shape(2)

    filename = 'BestLambdas'//f_name
    inquire(file=trim(filename), exist=exists)
    if (exists) then
        open(newunit = unit, file=trim(filename), action="read", iostat=ierr)
        if (ierr .eq. 0) then
            do i_waves=1,n_waves
                read(unit, *, iostat=readerr) param(:,i_waves)
                if (readerr .ne. 0) then
                    print*, "Error ", readerr, & 
                            " attempting to read file ", & 
                            trim(filename), "In wave number", i_waves
                    stop
                endif
                read(unit, *, iostat=readerr) 
            enddo
        else
            print*, "Error ", ierr ," attempting to open file ", trim(filename)
            stop
        endif
        close(unit)
    else
        print*, "Error -- cannot find file: ", trim(filename)
        stop
    endif

end subroutine read_parameters

subroutine write_momentum_dependence(param_samples,dr,dr_tail,r_max)
    implicit none
    real(dp), intent(in) :: param_samples(:,:,:), dr, dr_tail,r_max
    integer :: unit_st,ierr_st, unit_poly,ierr_poly
    character(len=128) :: filename, filename_poly
    real(dp) :: momentum, V_mean_st(1:n_q_operators),V_variance_st(1:n_q_operators)
    real(dp) :: V_mean_poly(1:n_av18q_operators),V_variance_poly(1:n_av18q_operators)
    integer :: i
    
    filename = 'momentum_dependence_st'//f_name
    filename_poly = 'momentum_dependence_poly'//f_name

    open(newunit = unit_st, file=trim(filename), action="write", iostat=ierr_st)
    open(newunit = unit_poly, file=trim(filename_poly), action="write", iostat=ierr_poly)
    if (ierr_st .eq. 0 .and. ierr_poly .eq. 0) then
        momentum = 1.e-1_dp
        do
            call  mc_momentum_dependence(param_samples,momentum,dr,dr_tail,r_max,V_mean_poly,V_mean_st,&
                & V_variance_poly,V_variance_st)
            write(unit_st,'(49e21.8)') momentum, (V_mean_st(i), V_variance_st(i),i=1,n_q_operators)
            write(unit_poly,'(49e21.8)') momentum, (V_mean_poly(i), V_variance_poly(i),i=1,n_av18q_operators)
            momentum = momentum + 1.e-1_dp
            if (momentum.gt.15._dp) exit
        enddo
    else
        if(ierr_st.ne.0) then
            print*, "Error ", ierr_st ," attempting to open file ", trim(filename)
        endif
        if(ierr_poly.ne.0) then
            print*, "Error ", ierr_poly ," attempting to open file ", trim(filename_poly)
        endif
        stop
    endif
    close(unit_st)
    close(unit_poly)
end subroutine write_momentum_dependence

subroutine write_av18_momentum(delta_r,r_max)
    implicit none
    real(dp), intent(in) :: delta_r,r_max
    real(dp), allocatable :: radii(:), av18_lambdas(:,:)
    character(len=128) :: filename, filename_st
    integer :: unit, ierr, unit_st,ierr_st
    real(dp) :: momentum
    real(dp) :: Vq(1:n_q_operators), V_poly(1:n_av18q_operators)
    integer :: i

    call sample_av18(delta_r,r_max,radii,av18_lambdas)
    
    filename = 'av18_momentum_dependence_poly_basis.dat'
    filename_st = 'av18_momentum_dependence_st_basis.dat'
    open(newunit = unit, file=trim(filename), action="write", iostat=ierr)
    open(newunit = unit_st, file=trim(filename_st), action="write", iostat=ierr_st)
    if (ierr .eq. 0 .and. ierr_st .eq. 0) then
        momentum = 1.e-1_dp
        do
            call transform_all_oper(momentum,av18_lambdas,radii,V_poly,Vq)
            write(unit,'(25e21.8)') momentum, V_poly
            write(unit_st,'(25e21.8)') momentum, Vq
            momentum = momentum + 1.e-1_dp
            if (momentum.gt.15._dp) exit
        enddo
    else
        if(ierr .ne. 0) then
            print*, "Error ", ierr ," attempting to open file ", trim(filename)
        endif 
        if(ierr_st .ne. 0) then
            print*, "Error ", ierr_st ," attempting to open file ", trim(filename_st)
        endif 
        stop
    endif
    close(unit)
    close(unit_st)

    
end subroutine write_av18_momentum

end module c2m_io