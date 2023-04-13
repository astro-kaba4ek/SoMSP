! Created by odawing on 11.09.22.

module scattering_calculation
    use regime
    use utils
    use contexts
    use constants
    use spheroidal_tmatrix
    use spheroidal_scattering
    use tmatrix_conversion
    use mode_functors
    use spheroidal_indicatrix
    use calculation_models
    use logging
    implicit none

    private
    public :: calculate_indicatrix, calculate_m

contains

    function calculate_indicatrix(scattering_context, minm, maxm, model, lnum, spherical_lnum, scatmatr_file) result(res)
        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: minm, maxm, lnum
        integer, optional, intent(in) :: spherical_lnum
        character(*), intent(in) :: model

        type(Node), allocatable :: queue(:)
        type(ModeCalculationResult), allocatable :: mode_res(:)
        type(ScatteringResult) :: res
        
        real(knd) :: accuracy, dtheta, dphi, value, bucket(4,4)
        integer :: qlen, i, m, sph_lnum, real_maxm, j, dim
        complex(knd), dimension(2 * lnum, minm:maxm) :: solution_te, solution_tm
        integer :: basis_types(minm:maxm)
        integer :: ntheta, nphi
        complex(knd), dimension(2,2,scattering_context%directions%nphi,scattering_context%directions%ntheta) :: ampl, update
        logical :: need_indicatrix

        real(knd), parameter :: low = 1.0e-8_knd
        type(SpheroidalCalculation) :: direction_calculation
        class(CalculationModel), allocatable :: typed_model
        type(SolutionForIndicatrix), allocatable :: solutions(:)
        character(*), intent(in) :: scatmatr_file

        integer :: theta_bucket_size, theta_bucket_start, theta_bucket_end, current_size, current_end

        call res%initialize()

        99 format('#',1A5,' ', 1A12,' ',6A24)
            write(*,99) 'm', 'potentials', &
            'Q_{TM}^{ext}', 'Q_{TM}^{sca}', 'Q_{TM}^{abs}', &
            'Q_{TE}^{ext}', 'Q_{TE}^{sca}', 'Q_{TE}^{abs}'

        sph_lnum = lnum
        if (present(spherical_lnum)) then
            sph_lnum = spherical_lnum
        endif

        ntheta = scattering_context%directions%ntheta
        nphi = scattering_context%directions%nphi

        ! model = 'all_uv'

        need_indicatrix = nphi > 0 .and. ntheta > 0 .and. &
        (model == 'all_uv' .or. model == 'uv_pq' .or. model == 'uv_pq_te_from_tm')

        if (need_indicatrix) open(SCAT_MATR_FD, file=trim(scatmatr_file), status='replace')

        ampl = 0

        if (nphi == 0) then
            theta_bucket_size = 0
        else
            theta_bucket_size = min(ntheta, min(GLOBAL_BUCKET / lnum, GLOBAL_BUCKET / nphi))
        endif

        theta_bucket_start = 1
        theta_bucket_end = theta_bucket_size

        typed_model = CalculationModel(model)
        solution_tm = 0
        solution_te = 0
        if (LOG_INFO) write(LOG_FD,*) '{INFO} start calculation with model '//model
        do m = minm, maxm
            call typed_model%build_mode_queue(m, lnum, spherical_lnum, queue)
            if (LOG_INFO) write(LOG_FD,*) 
            qlen = size(queue)
            if (qlen == 0) then
                cycle
            endif
            call log_node_queue(queue)

            accuracy = 0
            mode_res = calculate_m(scattering_context, queue, m, lnum, sph_lnum)

            do i = 1, qlen
                if (queue(i)%to_res) then
                    accuracy = max(accuracy, res%update_and_get_accuracy(queue(i)%info, mode_res(i)%factors))
                endif
            enddo

            call typed_model%print_mode_row(m, mode_res)

            if (need_indicatrix) then
                solutions = typed_model%solution_places(m)
                do i = 1, size(solutions)
                    dim = size(mode_res(solutions(i)%source_tm)%solution)
                    solution_tm(:dim,solutions(i)%m) = mode_res(solutions(i)%source_tm)%solution
                    solution_te(:dim,solutions(i)%m) = mode_res(solutions(i)%source_te)%solution
                    basis_types(solutions(i)%m) = solutions(i)%basis_type

                update = calculate_amplitude_matrix_m(solutions(i)%basis_type, scattering_context, solutions(i)%m, lnum, &
                theta_bucket_size, scattering_context%directions%thetas(theta_bucket_start:theta_bucket_end), &
                nphi, scattering_context%directions%phis, &
                direction_calculation, mode_res(solutions(i)%source_te)%solution, mode_res(solutions(i)%source_tm)%solution)
                accuracy = max(accuracy, update_amplitude_and_get_accuracy(ampl, update, ntheta, nphi))
                enddo

            endif

            if (size(queue) > 0 .and. accuracy < MIN_M_RATIO) then
                real_maxm = m
                exit
            endif
        end do
        deallocate(queue, mode_res)

        if (.not. need_indicatrix) then
            return
        endif


        100 format('#',2A8,' ',6A24)
        101 format(' ',2F8.2,' ',6F24.15)
        write(SCAT_MATR_FD,100) 'theta', 'phi', &
        'F_{11}', 'F_{21}', 'F_{33}', 'F_{43}'
        call print_scattering_matrix_bucket(&
        theta_bucket_size, scattering_context%directions%thetas(theta_bucket_start:theta_bucket_end), &
        nphi, scattering_context%directions%phis, ampl)
        do theta_bucket_start = theta_bucket_end + 1, ntheta, theta_bucket_size
            current_end = min(theta_bucket_start + theta_bucket_size - 1, ntheta)
            current_size = current_end - theta_bucket_start + 1
            ampl = 0
            do m = minm, maxm
                ampl(:,:,:,1:current_end) = ampl(:,:,:,1:current_end) + &
                calculate_amplitude_matrix_m(basis_types(m), scattering_context, m, lnum, &
                current_size, scattering_context%directions%thetas(theta_bucket_start:current_end), &
                nphi, scattering_context%directions%phis, &
                direction_calculation, solution_te(:,m), solution_tm(:,m))
            enddo

            call print_scattering_matrix_bucket(&
        current_size, scattering_context%directions%thetas(theta_bucket_start:current_end), &
        nphi, scattering_context%directions%phis, ampl(:,:,:,1:current_end))
        enddo

        deallocate(typed_model)
        if (allocated(solutions)) deallocate(solutions)
        if (need_indicatrix) close(SCAT_MATR_FD)
    end function calculate_indicatrix

    subroutine print_scattering_matrix_bucket(ntheta, thetas, nphi, phis, ampl)
        integer, intent(in) :: ntheta, nphi
        type(AngleType), intent(in) :: thetas(ntheta), phis(nphi)
        complex(knd), dimension(2,2,nphi, ntheta) :: ampl
        real(knd) :: bucket(4,4), start, finish
        integer :: i, j
        101 format(' ',2F8.2,' ',4F24.15)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate scattering matrix for bucket'
        call cpu_time(start)
        do i = 1, ntheta
            do j = 1, nphi
                bucket = get_scattering_matrix(ampl(:,:,j,i))
                write(SCAT_MATR_FD, 101) thetas(i)%value * 180q0 / PI, &
                    phis(j)%value * 180q0 / PI, &
                        bucket(1,1), bucket(2,1), bucket(3,3), bucket(4,3)
                ! write(*,'(4F24.15)') thetas(i)%value * 180q0 / PI, &
                ! phis(j)%value * 180q0 / PI, bucket
            end do
        end do
        call cpu_time(finish)
        call log_time('scattering matrix bucket', finish - start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate scattering matrix for bucket'
    end subroutine print_scattering_matrix_bucket

    function calculate_m(scattering_context, queue, m, lnum, spherical_lnum) result(results)
        type(ScatteringContext), intent(in) :: scattering_context
        type(Node), intent(in) :: queue(:)
        type(ModeCalculationResult), allocatable :: results(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        type(ComputationContext) :: computation_context

        integer :: len, i

        len = size(queue)

        if (allocated(results) .and. (size(results) /= len)) then
            deallocate(results)
        endif
        if (.not. allocated(results)) then
            allocate(results(len))
        endif

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate spheroidal functions'
        call computation_context%initialize(m, lnum, spherical_lnum, scattering_context)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate spheroidal functions'

        call calculate_mode_chain(scattering_context, computation_context, queue, results)

    end function calculate_m

    subroutine calculate_mode_chain(scattering_context, computation_context, queue, results)
        type(ScatteringContext), intent(in) :: scattering_context
        type(ComputationContext), intent(inout) :: computation_context
        type(Node), intent(in) :: queue(:)
        type(ModeCalculationResult) :: results(:)

        integer :: i

        do i = 1, size(queue)
            if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate node ', i
            if (size(queue(i)%previous) == 0) then
                call calculate_mode(queue(i), scattering_context, computation_context, results(i))
            elseif (size(queue(i)%previous) == 1) then
                call calculate_mode(queue(i), scattering_context, computation_context, results(i), &
                queue(queue(i)%previous(1)), results(queue(i)%previous(1))%tmatrix)
            else
                call assert(queue(i)%info == MODE_FAR_TETM .and. size(queue(i)%previous) == 2 .and. &
                    queue(queue(i)%previous(1))%info == MODE_FAR_TE_PQ .and.  &
                    queue(queue(i)%previous(2))%info == MODE_FAR_TM_PQ, &
                    'only far tm and te modes are united into one to transform to barber')
                call combine_tmatrices(results(queue(i)%previous(1))%tmatrix, results(queue(i)%previous(2))%tmatrix, &
                    results(i)%tmatrix)
            endif

            call log_matrix(queue(i)%info%to_string(), results(i)%tmatrix)
            if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate node ', i
        enddo
    end subroutine calculate_mode_chain

    subroutine combine_tmatrices(first_tmatrix, second_tmatrix, total_tmatrix)
        complex(knd), intent(in) :: first_tmatrix(:,:), second_tmatrix(:,:)
        complex(knd), allocatable, intent(out) :: total_tmatrix(:,:)

        integer :: first_size, second_size, total_size

        call assert(size(first_tmatrix(:,1)) == size(first_tmatrix(1,:)), 'first tmatrix is not square')
        call assert(size(second_tmatrix(:,1)) == size(second_tmatrix(1,:)), 'second tmatrix is not square')

        first_size = size(first_tmatrix(:,1))
        second_size = size(second_tmatrix(:,1))
        total_size = first_size + second_size

        if (allocated(total_tmatrix)) deallocate(total_tmatrix)
        allocate(total_tmatrix(total_size, total_size))

        total_tmatrix = 0
        total_tmatrix(1:first_size, 1:first_size) = first_tmatrix
        total_tmatrix(first_size + 1:total_size, first_size + 1: total_size) = second_tmatrix
    end subroutine combine_tmatrices

    subroutine calculate_mode(current_node, global_context, base_context, res, &
        old_node, old_tmatrix)
        type(Node), intent(in) :: current_node
        type(ScatteringContext), intent(in) :: global_context
        type(ComputationContext), intent(inout) :: base_context
        type(Node), optional, intent(in) :: old_node
        complex(knd), optional, intent(in) :: old_tmatrix(:,:)
        type(ModeCalculationResult), intent(out) :: res

        complex(knd), allocatable :: initial(:)
        type(ModeFunctors) :: funcs

        call res%initialize(current_node%get_matrix_size())

        call calculate_mode_tmatrix(global_context, base_context, current_node, res%tmatrix, old_node, old_tmatrix)


        if (current_node%need_calc) then
            funcs = get_mode_functions(current_node%info)
            ! allocates initial
            call funcs%initial_calculator(base_context, current_node%item, initial)
            res%solution = matmul(res%tmatrix, initial)
! 
            res%factors%Qext = funcs%extinction_calculator(base_context, current_node%item, res%solution)
            res%factors%Qsca = funcs%scattering_calculator(base_context, current_node%item, res%solution)
            
            deallocate(initial)
        endif
        
    end subroutine calculate_mode

    subroutine calculate_tmatrix_directly(scattering_context, computation_context, new_mode, new_tmatrix)
        type(ScatteringContext), intent(in) :: scattering_context
        type(ComputationContext), intent(inout) :: computation_context
        type(Node), intent(in) :: new_mode
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(SpheroidalContext), pointer :: context

        if (new_mode%info%basis /= SPHEROIDAL_BASIS) then
            call assert(.false., 'can only calculate sph directly')
        endif

        if (LOG_INFO) write(LOG_FD, *) '{INFO} calculate tmatrix directly '//trim(new_mode%info%to_string())&
        //' {'//trim(new_mode%item%to_string())//'}'
        if (new_mode%info%basis_type == UV) then
            context => computation_context%get_spheroidal_context(3, new_mode%item%m, new_mode%item%lnum, new_mode%item%lnum)
            call calculate_tmatrix_spheroidal_uv(&
                scattering_context, &
                context, &
                new_mode%info == MODE_SPH_TE_UV, &
                new_tmatrix)
        else
            context => computation_context%get_spheroidal_context(2, new_mode%item%m, new_mode%item%lnum, new_mode%item%lnum)
            call calculate_tmatrix_sph_pq(&
                scattering_context, &
                context, &
                new_mode%info == MODE_SPH_TE_PQ, &
                new_tmatrix)
        endif

    end subroutine calculate_tmatrix_directly

    subroutine calculate_mode_tmatrix(scattering_context, computation_context, new_node, new_tmatrix, old_node, old_tmatrix)
        type(ScatteringContext), intent(in) :: scattering_context
        type(ComputationContext), intent(inout) :: computation_context
        type(Node), intent(in) :: new_node
        complex(knd), intent(out) :: new_tmatrix(:,:)
        type(Node), optional, intent(in) :: old_node
        complex(knd), optional, intent(in) :: old_tmatrix(:,:)

        call assert(present(old_node) .eqv. present(old_tmatrix), 'only one of old mode and tmatrix is present')

        if (.not. present(old_node)) then
            call calculate_tmatrix_directly(scattering_context, computation_context, new_node, new_tmatrix)
        else
            call convert_tmatrix(computation_context, old_node, new_node, old_tmatrix, new_tmatrix)
        endif

    end subroutine calculate_mode_tmatrix

    
end module scattering_calculation