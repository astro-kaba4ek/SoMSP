program factors
    use regime
    use scattering_calculation
    use utils
    use constants
    use contexts
    use spheroidal_scatterer
    use spheroidal_indicatrix
    implicit none

    integer :: f, nol, matrix_size, spherical_lnum, minm, maxm, ntheta, nphi
    real(knd), allocatable :: rv(:), xv(:), ab(:)
    real(knd) :: alpha, lambda, theta0, theta1, phi0, phi1
    complex(knd), allocatable :: ri(:)
    ! type(ScatteringQuery) :: query
    type(ScatteringResult):: result
    type(ScatteringContext) :: global_context
    type(SpheroidalShape) :: shape
    character(32) :: model
    character(1024) :: input_file, scatmatr_file, arg
    integer :: i,j,k
    logical :: need_far

    if (LOG_INFO) open(LOG_FD, FILE=LOG_FILENAME, status='replace')
    input_file = 'input.txt'
    scatmatr_file = 'scattering_matrix.txt'
    do i = 1, command_argument_count(), 2
        call assert(i < command_argument_count(), 'Number of arguments must be even')
        call get_command_argument(i, arg)
        if (trim(arg) == '--input') then
            call get_command_argument(i + 1, input_file) 
        elseif (trim(arg) == '--scat-matr') then
            call get_command_argument(i + 1, scatmatr_file) 
        else
            call assert(.false., 'unknown parameter '//trim(arg))
        endif
    enddo

    call read_input(input_file, f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, model, &
    ntheta, theta0, theta1, nphi, phi0, phi1)

    call global_context%initialize(f, nol, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, &
    ntheta, theta0, theta1, nphi, phi0, phi1)

    result = calculate_indicatrix(global_context, minm, maxm, model, matrix_size, spherical_lnum, scatmatr_file)

    call shape%set(f, rv(1), ab(1), alpha)

    need_far = (model == 'uv_pq_te_from_tm_with_far')
    call log_mode_factors('Q SPH_TM', result%sph_tm)
    if (need_far) call log_mode_factors('Q FAR_TM', result%far_tm)
    call log_mode_factors('Q SPH_TE', result%sph_te)
    if (need_far) call log_mode_factors('Q FAR_TE', result%far_te)
    call log_mode_factors('C_SPH_TM', get_c_factors_from_q(result%sph_tm, shape))
    if (need_far) call log_mode_factors('C FAR_TM', get_c_factors_from_q(result%far_tm, shape))
    call log_mode_factors('C  SPH_TE', get_c_factors_from_q(result%sph_te, shape))
    if (need_far) call log_mode_factors('C FAR_TE', get_c_factors_from_q(result%far_te, shape))
    call log_mode_factors('C_norm SPH_TM', get_normalized_c_factors_from_q(result%sph_tm, shape))
    if (need_far) call log_mode_factors('C_norm FAR_TM', get_normalized_c_factors_from_q(result%far_tm, shape))
    call log_mode_factors('C norm SPH_TE', get_normalized_c_factors_from_q(result%sph_te, shape))
    if (need_far) call log_mode_factors('C norm FAR_TE', get_normalized_c_factors_from_q(result%far_te, shape))

    deallocate(rv, xv, ab, ri)
    if (LOG_INFO) close(LOG_FD)
end program factors
