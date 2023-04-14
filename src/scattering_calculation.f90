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

        ! MPI variables
        logical, allocatable :: node_logic(:)
        logical :: skip(Size_mpi-1), accuracy_logic
        integer, allocatable :: node_int(:), node_prev(:), nodes_previous(:,:), ij(:,:), MCR_ind(:), qlen_m(:)
        integer :: sum_size, sum_size2, k, s, c, c1, m_array(Size_mpi-1), iter
        real(knd), allocatable :: MCR_real(:)
        complex(knd), allocatable :: MCR_comp_mat(:), MCR_comp_arr(:)
        type(Node), allocatable :: queue_m(:,:)
        type(ModeCalculationResult), allocatable :: mode_res_m(:,:)


        call res%initialize()

        if (Rank == 0) then
            99 format('#',1A5,' ', 1A12,' ',6A24)
                write(*,99) 'm', 'potentials', &
                'Q_{TM}^{ext}', 'Q_{TM}^{sca}', 'Q_{TM}^{abs}', &
                'Q_{TE}^{ext}', 'Q_{TE}^{sca}', 'Q_{TE}^{abs}'
        endif
        
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


        accuracy_logic = .false.

        do iter = minm, maxm, Size_mpi-1
            
            m_array = [(m, m=iter, iter+Size_mpi-1-1)]

            if (Rank /= 0) then
                m = m_array(Rank)

                call typed_model%build_mode_queue(m, lnum, spherical_lnum, queue)
                if (LOG_INFO) write(LOG_FD,*) 

                qlen = size(queue)
                if (qlen == 0 .or. m > maxm) then

                    qlen = 0

                    call MPI_Ssend(qlen, 1, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)
                    call MPI_Barrier(MPI_COMM_WORLD, Err)
        
                else

                    call log_node_queue(queue)

                    mode_res = calculate_m(scattering_context, queue, m, lnum, sph_lnum)
                    
                    ! convert "queue" into 3 one-dimensional arrays and pass them
                    sum_size = qlen
                    do i = 1, qlen
                        sum_size = sum_size + size(queue(i)%previous)
                    enddo

                    allocate(node_logic(2*qlen), node_int((2+4)*qlen), node_prev(sum_size))

                    node_logic = [(queue(i)%need_calc, i = 1, qlen), (queue(i)%to_res, i = 1, qlen)]
                    node_int = [integer::]
                    node_prev = [integer::]
                    do i = 1, qlen
                        node_int = [node_int, queue(i)%info%basis, queue(i)%info%tmode, queue(i)%info%basis_type, queue(i)%info%num]
                        node_int = [node_int, queue(i)%item%m, queue(i)%item%lnum]
                        node_prev = [node_prev, size(queue(i)%previous), queue(i)%previous]
                    enddo

                    call MPI_Ssend(qlen, 1, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)
                    call MPI_Barrier(MPI_COMM_WORLD, Err)

                    call MPI_Ssend(node_logic, 2*qlen, MPI_LOGICAL, 0, m, MPI_COMM_WORLD, Err)
                    call MPI_Ssend(node_int, (2+4)*qlen, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)
                    call MPI_Ssend(node_prev, sum_size, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)

                    deallocate(node_logic, node_int, node_prev, queue)

                    
                    allocate(ij(qlen,2))
                    do i = 1, qlen
                        ij(i,1) = size(mode_res(i)%tmatrix, dim = 1) ! number of rows = length of columns
                        ij(i,2) = size(mode_res(i)%tmatrix, dim = 2) ! number of columns = length of rows
                    enddo

                    s = 0 ! number of elements in all matrices
                    sum_size2 = 0
                    do i = 1, qlen
                        s = s + product(ij(i,:))
                        sum_size2 = sum_size2 + size(mode_res(i)%solution)
                    enddo

                    allocate(MCR_real(2*qlen), MCR_comp_mat(s), MCR_comp_arr(sum_size2), MCR_ind(3*qlen))

                    MCR_real = [real(knd)::]
                    MCR_comp_mat = [complex(knd)::]
                    MCR_comp_arr = [complex(knd)::]
                    MCR_ind = [integer::]

                    do i = 1, qlen
                        MCR_real = [MCR_real, mode_res(i)%factors%Qext, mode_res(i)%factors%Qsca]
                        MCR_comp_mat = [MCR_comp_mat, (mode_res(i)%tmatrix(:,j), j=1, ij(i,1))]
                        MCR_comp_arr = [MCR_comp_arr, mode_res(i)%solution]
                        MCR_ind = [MCR_ind, ij(i,:), size(mode_res(i)%solution)]
                    enddo

                    if (knd == 16) then
                        call MPI_Ssend(MCR_real, 2*qlen, MPI_REAL16, 0, m, MPI_COMM_WORLD, Err)
                        call MPI_Ssend(MCR_comp_mat, s, MPI_COMPLEX16, 0, m, MPI_COMM_WORLD, Err)
                        call MPI_Ssend(MCR_comp_arr, sum_size2, MPI_COMPLEX16, 0, m, MPI_COMM_WORLD, Err)
                    elseif (knd == 8) then
                        call MPI_Ssend(MCR_real, 2*qlen, MPI_REAL8, 0, m, MPI_COMM_WORLD, Err)
                        call MPI_Ssend(MCR_comp_mat, s, MPI_COMPLEX8, 0, m, MPI_COMM_WORLD, Err)
                        call MPI_Ssend(MCR_comp_arr, sum_size2, MPI_COMPLEX8, 0, m, MPI_COMM_WORLD, Err)
                    endif
                    call MPI_Ssend(MCR_ind, 3*qlen, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)

                    deallocate(ij, MCR_real, MCR_comp_mat, MCR_comp_arr, MCR_ind, mode_res)

                endif

                ! checking the achieved accuracy and exiting the outer loop
                ! call MPI_Barrier(MPI_COMM_WORLD, Err)
                call MPI_Bcast(accuracy_logic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Err)
                if (accuracy_logic) exit 


            else!if (Rank == 0) then

                ! collect an array with information about the lengths of the "queue" and "mode_res" arrays 
                ! (they may differ at each iteration of the loop)
                skip = .false.
                allocate(qlen_m(iter:iter+Size_mpi-1-1)); qlen_m = 0
                do k = 1, size(m_array)

                    call MPI_Recv(qlen, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, Status, Err)

                    ! get the thread tag, that is, "m"
                    m = Status%MPI_TAG
                    qlen_m(m) = qlen

                    if (qlen == 0) skip(k) = .true. ! TRUE == skip this "m"

                enddo
                call MPI_Barrier(MPI_COMM_WORLD, Err)


                allocate(mode_res_m(maxval(qlen_m),iter:iter+Size_mpi-1-1), & 
                        queue_m(maxval(qlen_m),iter:iter+Size_mpi-1-1))

                do k = 1, size(m_array)-count(skip)
                    
                    ! get "m"
                    call MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, Status, Err) 
                    m = Status%MPI_TAG
                    qlen = qlen_m(m)

                    ! get encrypted "queue"
                    allocate(node_logic(2*qlen))
                    call MPI_Recv(node_logic, 2*qlen, MPI_LOGICAL, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                    allocate(node_int((2+4)*qlen))
                    call MPI_Recv(node_int, (2+4)*qlen, MPI_INTEGER, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                    call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                    call MPI_Get_count(Status, MPI_INTEGER, sum_size, Err)

                    allocate(node_prev(sum_size))
                    call MPI_Recv(node_prev, sum_size, MPI_INTEGER, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                    ! decoding
                    c = 1; c1 = 0
                    allocate(queue(qlen))
                    do i = 1, qlen
                        queue(i)%need_calc = node_logic(i)
                        queue(i)%to_res = node_logic(i+qlen)

                        queue(i)%info%basis = node_int(6*(i-1)+1)
                        queue(i)%info%tmode = node_int(6*(i-1)+2)
                        queue(i)%info%basis_type = node_int(6*(i-1)+3)
                        queue(i)%info%num = node_int(6*(i-1)+4)

                        queue(i)%item%m = node_int(6*(i-1)+5)
                        queue(i)%item%lnum = node_int(6*(i-1)+6)

                        c1 = node_prev(c)
                        queue(i)%previous = node_prev(c+1:c+c1)
                        c = c + c1 + 1
                    enddo
                    
                    queue_m(:qlen,m) = queue

                    deallocate(node_prev, node_logic, node_int, queue)


                    if (knd == 16) then
                        allocate(MCR_real(2*qlen))
                        call MPI_Recv(MCR_real, 2*qlen, MPI_REAL16, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
            
                        call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                        call MPI_Get_count(Status, MPI_COMPLEX16, s, Err)

                        allocate(MCR_comp_mat(s))
                        call MPI_Recv(MCR_comp_mat, s, MPI_REAL16, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
            
                        call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                        call MPI_Get_count(Status, MPI_COMPLEX16, sum_size2, Err)

                        allocate(MCR_comp_arr(sum_size2))
                        call MPI_Recv(MCR_comp_arr, sum_size2, MPI_REAL16, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                    elseif (knd == 8) then
                        allocate(MCR_real(2*qlen))
                        call MPI_Recv(MCR_real, 2*qlen, MPI_REAL8, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
            
                        call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                        call MPI_Get_count(Status, MPI_COMPLEX8, s, Err)

                        allocate(MCR_comp_mat(s))
                        call MPI_Recv(MCR_comp_mat, s, MPI_COMPLEX8, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
            
                        call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                        call MPI_Get_count(Status, MPI_COMPLEX8, sum_size2, Err)

                        allocate(MCR_comp_arr(sum_size2))
                        call MPI_Recv(MCR_comp_arr, sum_size2, MPI_COMPLEX8, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                    endif

                    allocate(MCR_ind(3*qlen))
                    call MPI_Recv(MCR_ind, 3*qlen, MPI_INTEGER, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                    allocate(mode_res(qlen))
                    allocate(ij(qlen,2))
                    c1 = 0
                    do i = 1, qlen
                        mode_res(i)%factors%Qext = MCR_real(2*i-1)
                        mode_res(i)%factors%Qsca = MCR_real(2*i)

                        ij(i,:) = MCR_ind(3*i-2:3*i-1)
                        c = MCR_ind(3*i)

                        mode_res(i)%solution = MCR_comp_arr(i*c-c+1:i*c)

                        mode_res(i)%tmatrix = reshape(MCR_comp_mat(c1+1:c1+product(ij(i,:))), [ij(i,1), ij(i,2)])
                        c1 = c1 + product(ij(i,:))
                    enddo

                    mode_res_m(:qlen,m) = mode_res

                    deallocate(ij, MCR_real, MCR_comp_mat, MCR_comp_arr, MCR_ind, mode_res)

                enddo


                ! essentially the same loop as it was, 
                ! but all the heavy calculations have already been calculated and collected above
                do k = 1, size(m_array)
                    m = m_array(k)
                    qlen = qlen_m(m)

                    if (qlen == 0) cycle

                    accuracy = 0

                    do i = 1, qlen
                        if (queue_m(i,m)%to_res) then
                            accuracy = max(accuracy, res%update_and_get_accuracy(queue_m(i,m)%info, mode_res_m(i,m)%factors))
                        endif
                        ! print*, "kekeke",m, accuracy
                    enddo

                    call typed_model%print_mode_row(m, mode_res_m(:qlen,m))

                    if (need_indicatrix) then
                        solutions = typed_model%solution_places(m)
                        do i = 1, size(solutions)
                            dim = size(mode_res_m(solutions(i)%source_tm,m)%solution)
                            solution_tm(:dim,solutions(i)%m) = mode_res_m(solutions(i)%source_tm,m)%solution
                            solution_te(:dim,solutions(i)%m) = mode_res_m(solutions(i)%source_te,m)%solution
                            basis_types(solutions(i)%m) = solutions(i)%basis_type
        
                            update = calculate_amplitude_matrix_m(solutions(i)%basis_type, scattering_context, solutions(i)%m, & 
                            lnum, theta_bucket_size, scattering_context%directions%thetas(theta_bucket_start:theta_bucket_end), & 
                            nphi, scattering_context%directions%phis, direction_calculation, & 
                            mode_res_m(solutions(i)%source_te,m)%solution, mode_res_m(solutions(i)%source_tm,m)%solution)
                            accuracy = max(accuracy, update_amplitude_and_get_accuracy(ampl, update, ntheta, nphi))
                            
                            ! print*, "ggggg", m, ampl
                            ! print*, "ggggg", m, update
                            ! print*, "ggggg", m, ntheta
                            ! print*, "ggggg", m, nphi
                            ! print*, "ggggg", m, update_amplitude_and_get_accuracy(ampl, update, ntheta, nphi)
                            ! print*, "lollol", m, accuracy
                        enddo
        
                    endif

                    print*, "accuracy", size(queue_m(:qlen,m)), accuracy, MIN_M_RATIO

                    if (size(queue_m(:qlen,m)) > 0 .and. accuracy < MIN_M_RATIO) then
                        accuracy_logic = .true.
                        real_maxm = m
                        exit
                    endif

                enddo

                deallocate(mode_res_m, queue_m, qlen_m)

                ! inform all threads about the achieved accuracy and exit the outer loop
                ! call MPI_Barrier(MPI_COMM_WORLD, Err)
                call MPI_Bcast(accuracy_logic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Err)
                if (accuracy_logic) exit 

            endif

        enddo 

        print*, "real_maxm", real_maxm, rank

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