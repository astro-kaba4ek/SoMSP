program factors
    use regime
    use scattering_calculation
    use utils
    use constants
    use contexts
    use spheroidal_scatterer
    use spheroidal_indicatrix
    ! use calculation_models
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

    ! integer :: m
    ! class(CalculationModel), allocatable :: typed_model

    ! ! MPI variables
    ! logical :: accuracy_logic
    ! real(knd) :: time_t1, time_t2, time_t3, time_t4
    character(3) :: str_rank



    ! MPI on
    call MPI_Init(Err)
    call MPI_Comm_rank(MPI_COMM_WORLD, Rank, Err)
    call MPI_Comm_size(MPI_COMM_WORLD, Size_mpi, Err)

    write(str_rank, "(i3)") Rank
    if (LOG_INFO) open(LOG_FD, FILE=trim(str_rank)//"n_"//trim(LOG_FILENAME), status='replace')

    if (Rank == 0) then
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

        call MPI_Bcast(input_file, len(input_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, Err)
        call MPI_Bcast(scatmatr_file, len(scatmatr_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, Err)
    
    else 
        call MPI_Bcast(input_file, len(input_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, Err)
        call MPI_Bcast(scatmatr_file, len(scatmatr_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, Err)
    endif

    call read_input(input_file, f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, model, &
    ntheta, theta0, theta1, nphi, phi0, phi1)

    call global_context%initialize(f, nol, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, &
    ntheta, theta0, theta1, nphi, phi0, phi1)

    ! result = calculate_indicatrix(global_context, minm, maxm, model, matrix_size, spherical_lnum, scatmatr_file)
    result = calculate_indicatrix_2(global_context, minm, maxm, model, matrix_size, spherical_lnum, scatmatr_file)
   
    ! call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, knd, MPI_REAL_knd, Err)
    ! call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_COMPLEX, 2*knd, MPI_COMPLEX_knd, Err)
    
    ! typed_model = CalculationModel(model)

    ! accuracy_logic = .false.

    ! ! allocate(m_array_all(minm:maxm))
    ! ! m_array_all = [(m, m=minm, maxm)]

    ! if (Rank == 0) then
    !     call for_rank_0(global_context, minm, maxm, model, matrix_size, scatmatr_file,&
    !      typed_model, accuracy_logic, result)
    ! else 
    !     call for_not_rank_0(global_context, minm, maxm, matrix_size, spherical_lnum, &
    !     typed_model, accuracy_logic)
    ! endif

    ! deallocate(typed_model)

    if (Rank == 0) then

        call shape%set(f, rv(1), ab(1), alpha)

        need_far = (model == 'uv_pq_te_from_tm_with_far')
        call log_mode_factors('Q_SPH_TM', result%sph_tm)
        if (need_far) call log_mode_factors('Q_FAR_TM', result%far_tm)
        call log_mode_factors('Q_SPH_TE', result%sph_te)
        if (need_far) call log_mode_factors('Q_FAR_TE', result%far_te)
        call log_mode_factors('C_SPH_TM', get_c_factors_from_q(result%sph_tm, shape))
        if (need_far) call log_mode_factors('C_FAR_TM', get_c_factors_from_q(result%far_tm, shape))
        call log_mode_factors('C_SPH_TE', get_c_factors_from_q(result%sph_te, shape))
        if (need_far) call log_mode_factors('C_FAR_TE', get_c_factors_from_q(result%far_te, shape))
        call log_mode_factors('C_norm_SPH_TM', get_normalized_c_factors_from_q(result%sph_tm, shape))
        if (need_far) call log_mode_factors('C_norm_FAR_TM', get_normalized_c_factors_from_q(result%far_tm, shape))
        call log_mode_factors('C_norm_SPH_TE', get_normalized_c_factors_from_q(result%sph_te, shape))
        if (need_far) call log_mode_factors('C_norm_FAR_TE', get_normalized_c_factors_from_q(result%far_te, shape))

        if (LOG_INFO) close(LOG_FD)
    end if

    deallocate(rv, xv, ab, ri)

    ! MPI off
    call MPI_Finalize(Err)

end program factors
