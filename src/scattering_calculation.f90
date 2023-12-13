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
    public :: calculate_indicatrix, calculate_m, calculate_indicatrix_2, for_not_rank_0, for_rank_0

contains

    subroutine for_not_rank_0(scattering_context, minm, maxm, lnum, spherical_lnum, typed_model, accuracy_logic)

        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: minm, maxm, lnum
        integer, optional, intent(in) :: spherical_lnum

        integer :: qlen, i, m, sph_lnum, j
        type(Node), allocatable :: queue(:)
        type(ModeCalculationResult), allocatable :: mode_res(:)
        class(CalculationModel), intent(in) :: typed_model

        ! MPI variables
        logical, intent(inout) :: accuracy_logic
        integer, allocatable :: ij(:,:)
        integer :: s, sum_size, m_array(Size_mpi-1), iter, m_array_all(minm:maxm), position, len_prev, queue_buf_len
        type(MPI_Datatype) :: MPI_ModeInfo_type, MPI_ModeItem_type, MPI_Node_type
        real(knd) :: time_t3, time_t4
        integer(MPI_ADDRESS_KIND) :: lb, MPI_ModeInfo_size, MPI_ModeItem_size
        character, allocatable :: queue_buf(:)


        sph_lnum = lnum
        if (present(spherical_lnum)) then
            sph_lnum = spherical_lnum
        endif
        call MPI_Barrier(MPI_COMM_WORLD, Err)

        m_array_all = [(m, m=minm, maxm)]


        do iter = minm, maxm, Size_mpi-1
            m_array = m_array_all(iter : iter+Size_mpi-1-1)
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

                call cpu_time(time_t3)

                mode_res = calculate_m(scattering_context, queue, m, lnum, sph_lnum)
            
                call cpu_time(time_t4)
                ! print*, rank, iter, "mode_res time", time_t4-time_t3
                
                sum_size = 0
                do i = 1, qlen
                    sum_size = sum_size + size(queue(i)%previous)
                enddo

                call MPI_Ssend(qlen, 1, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)
                call MPI_Barrier(MPI_COMM_WORLD, Err)

                queue_buf_len = 4*2*qlen + 4*4*qlen + 4*2*qlen + 4*1*qlen + 4*sum_size

                allocate(queue_buf(queue_buf_len))
                position = 0
                do i=1, qlen
                    len_prev = size(queue(i)%previous)
                    call MPI_Pack(queue(i)%need_calc, 1, MPI_LOGICAL, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%to_res, 1, MPI_LOGICAL, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%info%basis, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%info%tmode, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%info%basis_type, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%info%num, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%item%m, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%item%lnum, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(len_prev, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(queue(i)%previous, len_prev, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                enddo

                call MPI_Ssend(queue_buf, queue_buf_len, MPI_PACKED, 0, m, MPI_COMM_WORLD, Err)

                deallocate(queue_buf, queue)

                
                allocate(ij(2,qlen))
                do i = 1, qlen
                    ! ij(i,1) = size(mode_res(i)%tmatrix, dim = 1) ! number of rows = length of columns
                    ! ij(i,2) = size(mode_res(i)%tmatrix, dim = 2) ! number of columns = length of rows
                    ij(:,i) = shape(mode_res(i)%tmatrix)
                enddo

                s = sum(product(ij, dim=1)) ! number of elements in all matrices
                sum_size = 0
                do i = 1, qlen
                    sum_size = sum_size + size(mode_res(i)%solution)
                enddo

                queue_buf_len = knd*2*qlen + (2*knd)*s*qlen + (2*knd)*sum_size + 4*3*qlen

                allocate(queue_buf(queue_buf_len))
                position = 0
                do i=1, qlen
                    len_prev = size(mode_res(i)%solution)
                    call MPI_Pack(mode_res(i)%factors%Qext, 1, MPI_REAL_knd, queue_buf, queue_buf_len, &
                                    position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(mode_res(i)%factors%Qsca, 1, MPI_REAL_knd, queue_buf, queue_buf_len, &
                                    position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(ij(:,i), 2, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(len_prev, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(pack(mode_res(i)%tmatrix, .true.), product(ij(:,i)), MPI_COMPLEX_knd, queue_buf, queue_buf_len, &
                                    position, MPI_COMM_WORLD, Err)
                    call MPI_Pack(mode_res(i)%solution, len_prev, MPI_COMPLEX_knd, queue_buf, queue_buf_len, &
                                    position, MPI_COMM_WORLD, Err)
                enddo

                call MPI_Ssend(queue_buf, queue_buf_len, MPI_PACKED, 0, m, MPI_COMM_WORLD, Err)

                deallocate(queue_buf, ij, mode_res)
            endif

            ! checking the achieved accuracy and exiting the outer loop
            ! call MPI_Barrier(MPI_COMM_WORLD, Err)
            call MPI_Bcast(accuracy_logic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Err)
            if (accuracy_logic) exit

        enddo

    end subroutine for_not_rank_0

    subroutine for_rank_0(scattering_context, minm, maxm, model, lnum, scatmatr_file, typed_model, accuracy_logic, res)
        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: minm, maxm, lnum
        character(*), intent(in) :: model

        type(Node), allocatable :: queue(:)
        type(ModeCalculationResult), allocatable :: mode_res(:)
        type(ScatteringResult), intent(out) :: res

        real(knd) :: accuracy
        integer :: qlen, i, m, real_maxm, dim
        complex(knd), dimension(2 * lnum, minm:maxm) :: solution_te, solution_tm
        integer :: basis_types(minm:maxm)
        integer :: ntheta, nphi
        complex(knd), dimension(2,2,scattering_context%directions%nphi,scattering_context%directions%ntheta) :: ampl, update
        logical :: need_indicatrix

        type(SpheroidalCalculation) :: direction_calculation
        class(CalculationModel), intent(in) :: typed_model
        character(*), intent(in) :: scatmatr_file
        type(SolutionForIndicatrix), allocatable :: solutions(:)

        integer :: theta_bucket_size, theta_bucket_start, theta_bucket_end, current_size, current_end

        ! MPI variables
        logical :: skip(Size_mpi-1)
        logical, intent(inout) :: accuracy_logic
        integer, allocatable :: qlen_m(:), prevv(:)
        integer :: k, m_array(Size_mpi-1), iter, m_array_all(minm:maxm), ij_2(2), position, len_prev, queue_buf_len
        complex(knd), allocatable :: MCR_comp_mat(:), MCR_comp_arr(:)
        type(Node), allocatable :: queue_m(:,:)
        type(ModeCalculationResult), allocatable :: mode_res_m(:,:)

        real(knd) :: time_t3, time_t4
        character, allocatable :: queue_buf(:)



        call res%initialize()

        99 format('#',1A5,' ', 1A12,' ',6A24)
            write(*,99) 'm', 'potentials', &
            'Q_{TM}^{ext}', 'Q_{TM}^{sca}', 'Q_{TM}^{abs}', &
            'Q_{TE}^{ext}', 'Q_{TE}^{sca}', 'Q_{TE}^{abs}'


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

        solution_tm = 0
        solution_te = 0
        if (LOG_INFO) write(LOG_FD,*) '{INFO} start calculation with model '//model

        call MPI_Barrier(MPI_COMM_WORLD, Err)

        m_array_all = [(m, m=minm, maxm)]
        
        do iter = minm, maxm, Size_mpi-1
            m_array = m_array_all(iter : iter+Size_mpi-1-1)
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

                call MPI_Get_count(Status, MPI_PACKED, queue_buf_len, Err)

                allocate(queue_buf(queue_buf_len))
                call MPI_Recv(queue_buf, queue_buf_len, MPI_PACKED, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                allocate(queue(qlen))
                position = 0
                do i=1, qlen
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%need_calc, 1, MPI_LOGICAL, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%to_res, 1, MPI_LOGICAL, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%basis, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%tmode, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%basis_type, 1, &
                                    MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%num, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%item%m, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%item%lnum, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, len_prev, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)

                    allocate(prevv(len_prev))
                    call MPI_Unpack(queue_buf, queue_buf_len, position, prevv, len_prev, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    queue(i)%previous = prevv
                    deallocate(prevv)
                enddo

                queue_m(:qlen,m) = queue
                deallocate(queue_buf, queue)


                call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                call MPI_Get_count(Status, MPI_PACKED, queue_buf_len, Err)

                allocate(queue_buf(queue_buf_len))
                call MPI_Recv(queue_buf, queue_buf_len, MPI_PACKED, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err)

                allocate(mode_res(qlen))
                position = 0
                do i=1, qlen
                    call MPI_Unpack(queue_buf, queue_buf_len, position, mode_res(i)%factors%Qext, 1, &
                                    MPI_REAL_knd, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, mode_res(i)%factors%Qsca, 1, &
                                    MPI_REAL_knd, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, ij_2, 2, MPI_INTEGER, MPI_COMM_WORLD, Err)
                    call MPI_Unpack(queue_buf, queue_buf_len, position, len_prev, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)

                    allocate(MCR_comp_mat(product(ij_2)))
                    call MPI_Unpack(queue_buf, queue_buf_len, position, MCR_comp_mat, product(ij_2), &
                                    MPI_COMPLEX_knd, MPI_COMM_WORLD, Err)
                    mode_res(i)%tmatrix = reshape(MCR_comp_mat, ij_2)
                    deallocate(MCR_comp_mat)
                    
                    allocate(MCR_comp_arr(len_prev))
                    call MPI_Unpack(queue_buf, queue_buf_len, position, MCR_comp_arr, len_prev, &
                                    MPI_COMPLEX_knd, MPI_COMM_WORLD, Err)
                    mode_res(i)%solution = MCR_comp_arr
                    deallocate(MCR_comp_arr)
                enddo

                mode_res_m(:qlen,m) = mode_res
                deallocate(queue_buf, mode_res)

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

                ! print*, "accuracy", size(queue_m(:qlen,m)), accuracy, MIN_M_RATIO

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
        enddo
    
        ! print*, "real_maxm", real_maxm, rank

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

        if (allocated(solutions)) deallocate(solutions)
        if (need_indicatrix) close(SCAT_MATR_FD)
        
    end subroutine for_rank_0

    function calculate_indicatrix_2(scattering_context, minm, maxm, model, lnum, spherical_lnum, scatmatr_file) result(res)
        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: minm, maxm, lnum
        integer, optional, intent(in) :: spherical_lnum
        character(*), intent(in) :: model

        type(ScatteringResult) :: res
        
        integer :: m
        class(CalculationModel), allocatable :: typed_model
        character(*), intent(in) :: scatmatr_file

        ! MPI variables
        logical :: accuracy_logic
        real(knd) :: time_t3, time_t4
      

        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, knd, MPI_REAL_knd, Err)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_COMPLEX, 2*knd, MPI_COMPLEX_knd, Err)
        
        typed_model = CalculationModel(model)

        accuracy_logic = .false.

        if (Rank == 0) then
            call for_rank_0(scattering_context, minm, maxm, model, lnum, scatmatr_file, &
                             typed_model, accuracy_logic, res)
        else 
            call for_not_rank_0(scattering_context, minm, maxm, lnum, spherical_lnum, &
                             typed_model, accuracy_logic)
        endif

        deallocate(typed_model)

    end function calculate_indicatrix_2

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
        logical :: skip(Size_mpi-1), accuracy_logic
        integer, allocatable :: ij(:,:), qlen_m(:), prevv(:)
        integer :: sum_size, k, s, m_array(Size_mpi-1), iter, ij_2(2), position, len_prev, queue_buf_len
        complex(knd), allocatable :: MCR_comp_mat(:), MCR_comp_arr(:)
        type(Node), allocatable :: queue_m(:,:)
        type(ModeCalculationResult), allocatable :: mode_res_m(:,:)

        real(knd) :: time_t3, time_t4
        character, allocatable :: queue_buf(:)
           
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, knd, MPI_REAL_knd, Err)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_COMPLEX, 2*knd, MPI_COMPLEX_knd, Err)


        if (Rank == 0) then

            call res%initialize()

            99 format('#',1A5,' ', 1A12,' ',6A24)
                write(*,99) 'm', 'potentials', &
                'Q_{TM}^{ext}', 'Q_{TM}^{sca}', 'Q_{TM}^{abs}', &
                'Q_{TE}^{ext}', 'Q_{TE}^{sca}', 'Q_{TE}^{abs}'


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

            solution_tm = 0
            solution_te = 0
            if (LOG_INFO) write(LOG_FD,*) '{INFO} start calculation with model '//model

        else 
            sph_lnum = lnum
            if (present(spherical_lnum)) then
                sph_lnum = spherical_lnum
            endif

        endif
        call MPI_Barrier(MPI_COMM_WORLD, Err)
        

        typed_model = CalculationModel(model)
       
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

                    call cpu_time(time_t3)
    
                    mode_res = calculate_m(scattering_context, queue, m, lnum, sph_lnum)
                
                    call cpu_time(time_t4)
                    ! print*, rank, iter, "mode_res time", time_t4-time_t3
                    
                    sum_size = 0
                    do i = 1, qlen
                        sum_size = sum_size + size(queue(i)%previous)
                    enddo
    
                    call MPI_Ssend(qlen, 1, MPI_INTEGER, 0, m, MPI_COMM_WORLD, Err)
                    call MPI_Barrier(MPI_COMM_WORLD, Err)
    
                    queue_buf_len = 4*2*qlen + 4*4*qlen + 4*2*qlen + 4*1*qlen + 4*sum_size
    
                    allocate(queue_buf(queue_buf_len))
                    position = 0
                    do i=1, qlen
                        len_prev = size(queue(i)%previous)
                        call MPI_Pack(queue(i)%need_calc, 1, MPI_LOGICAL, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%to_res, 1, MPI_LOGICAL, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%info%basis, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%info%tmode, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%info%basis_type, 1, MPI_INTEGER, queue_buf, queue_buf_len, &
                                        position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%info%num, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%item%m, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%item%lnum, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(len_prev, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(queue(i)%previous, len_prev, MPI_INTEGER, queue_buf, queue_buf_len, &
                                        position, MPI_COMM_WORLD, Err)
                    enddo
    
                    call MPI_Ssend(queue_buf, queue_buf_len, MPI_PACKED, 0, m, MPI_COMM_WORLD, Err)
    
                    deallocate(queue_buf, queue)
    
                    
                    allocate(ij(2,qlen))
                    do i = 1, qlen
                        ! ij(i,1) = size(mode_res(i)%tmatrix, dim = 1) ! number of rows = length of columns
                        ! ij(i,2) = size(mode_res(i)%tmatrix, dim = 2) ! number of columns = length of rows
                        ij(:,i) = shape(mode_res(i)%tmatrix)
                    enddo
    
                    s = sum(product(ij, dim=1)) ! number of elements in all matrices
                    sum_size = 0
                    do i = 1, qlen
                        sum_size = sum_size + size(mode_res(i)%solution)
                    enddo
    
                    queue_buf_len = knd*2*qlen + (2*knd)*s*qlen + (2*knd)*sum_size + 4*3*qlen
    
                    allocate(queue_buf(queue_buf_len))
                    position = 0
                    do i=1, qlen
                        len_prev = size(mode_res(i)%solution)
                        call MPI_Pack(mode_res(i)%factors%Qext, 1, MPI_REAL_knd, queue_buf, queue_buf_len, &
                                        position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(mode_res(i)%factors%Qsca, 1, MPI_REAL_knd, queue_buf, queue_buf_len, &
                                        position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(ij(:,i), 2, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(len_prev, 1, MPI_INTEGER, queue_buf, queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(pack(mode_res(i)%tmatrix, .true.), product(ij(:,i)), MPI_COMPLEX_knd, queue_buf, &
                                            queue_buf_len, position, MPI_COMM_WORLD, Err)
                        call MPI_Pack(mode_res(i)%solution, len_prev, MPI_COMPLEX_knd, queue_buf, queue_buf_len, &
                                        position, MPI_COMM_WORLD, Err)
                    enddo
    
                    call MPI_Ssend(queue_buf, queue_buf_len, MPI_PACKED, 0, m, MPI_COMM_WORLD, Err)
    
                    deallocate(queue_buf, ij, mode_res)

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

                    call MPI_Get_count(Status, MPI_PACKED, queue_buf_len, Err)

                    allocate(queue_buf(queue_buf_len))
                    call MPI_Recv(queue_buf, queue_buf_len, MPI_PACKED, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 

                    allocate(queue(qlen))
                    position = 0
                    do i=1, qlen
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%need_calc, 1, MPI_LOGICAL, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%to_res, 1, MPI_LOGICAL, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%basis, 1, &
                                        MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%tmode, 1, &
                                        MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%basis_type, 1, &
                                        MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%info%num, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%item%m, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, queue(i)%item%lnum, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, len_prev, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)

                        allocate(prevv(len_prev))
                        call MPI_Unpack(queue_buf, queue_buf_len, position, prevv, len_prev, MPI_INTEGER, MPI_COMM_WORLD, Err)
                        queue(i)%previous = prevv
                        deallocate(prevv)
                    enddo

                    queue_m(:qlen,m) = queue
                    deallocate(queue_buf, queue)


                    call MPI_Probe(MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err) 
                    call MPI_Get_count(Status, MPI_PACKED, queue_buf_len, Err)

                    allocate(queue_buf(queue_buf_len))
                    call MPI_Recv(queue_buf, queue_buf_len, MPI_PACKED, MPI_ANY_SOURCE, m, MPI_COMM_WORLD, Status, Err)

                    allocate(mode_res(qlen))
                    position = 0
                    do i=1, qlen
                        call MPI_Unpack(queue_buf, queue_buf_len, position, mode_res(i)%factors%Qext, 1, &
                                        MPI_REAL_knd, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, mode_res(i)%factors%Qsca, 1, &
                                        MPI_REAL_knd, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, ij_2, 2, MPI_INTEGER, MPI_COMM_WORLD, Err)
                        call MPI_Unpack(queue_buf, queue_buf_len, position, len_prev, 1, MPI_INTEGER, MPI_COMM_WORLD, Err)

                        allocate(MCR_comp_mat(product(ij_2)))
                        call MPI_Unpack(queue_buf, queue_buf_len, position, MCR_comp_mat, product(ij_2), &
                                        MPI_COMPLEX_knd, MPI_COMM_WORLD, Err)
                        mode_res(i)%tmatrix = reshape(MCR_comp_mat, ij_2)
                        deallocate(MCR_comp_mat)
                        
                        allocate(MCR_comp_arr(len_prev))
                        call MPI_Unpack(queue_buf, queue_buf_len, position, MCR_comp_arr, len_prev, &
                                        MPI_COMPLEX_knd, MPI_COMM_WORLD, Err)
                        mode_res(i)%solution = MCR_comp_arr
                        deallocate(MCR_comp_arr)
                    enddo

                    mode_res_m(:qlen,m) = mode_res
                    deallocate(queue_buf, mode_res)

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

                    ! print*, "accuracy", size(queue_m(:qlen,m)), accuracy, MIN_M_RATIO

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

        if (Rank == 0) then 

            ! print*, "real_maxm", real_maxm, rank

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

            if (allocated(solutions)) deallocate(solutions)
            if (need_indicatrix) close(SCAT_MATR_FD)

        endif

        deallocate(typed_model)

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