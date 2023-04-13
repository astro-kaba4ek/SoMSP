! Created by odawing on 11.09.22.

module utils
    use regime
    use constants
    use geometry
    implicit none

    real(knd), parameter :: LOWEST_UPDATE = 1.0e-128_knd

    type :: ModeInfo
        integer :: basis, tmode, basis_type, num
    contains
        procedure :: to_string => mode_info_to_string
    end type ModeInfo

    interface operator(.in.)
        procedure :: mode_operator_in
    end interface

    type ModeItem
        integer :: m, lnum
    contains
        procedure :: to_string => print_mode_item
    end type ModeItem
    ! used in arrays as a chain mode calculation request
    type Node
        ! which mode to calculate
        type(ModeInfo) :: info
        type(ModeItem) :: item
        ! number of the previous mode in the chain
        ! if this mode needs to be calculated directly, it is equal to the number of the current mode
        integer, allocatable :: previous(:)
        ! need to calculate the factors for the current mode
        logical :: need_calc
        ! used in final calculation
        logical :: to_res
    contains
        procedure :: get_matrix_size => get_matrix_size_by_node
        final :: delete_node
    end type Node

    interface Node
        procedure :: construct_node
    end interface Node

    ! type ModeTransition
    !     type(ModeInfo) :: old_mode, new_mode
    ! end type ModeTransition

    type(ModeInfo), parameter :: MODE_SPH_TE_UV = ModeInfo(basis=SPHEROIDAL_BASIS, tmode=TE, basis_type=UV, num=1)
    type(ModeInfo), parameter :: MODE_SPH_TM_UV = ModeInfo(basis=SPHEROIDAL_BASIS, tmode=TM, basis_type=UV, num=2)
    type(ModeInfo), parameter :: MODE_SPH_TE_PQ = ModeInfo(basis=SPHEROIDAL_BASIS, tmode=TE, basis_type=PQ, num=3)
    type(ModeInfo), parameter :: MODE_SPH_TM_PQ = ModeInfo(basis=SPHEROIDAL_BASIS, tmode=TM, basis_type=PQ, num=4)
    type(ModeInfo), parameter :: MODE_FAR_TE_UV = ModeInfo(basis=FAR_BASIS, tmode=TE, basis_type=UV, num=5)
    type(ModeInfo), parameter :: MODE_FAR_TM_UV = ModeInfo(basis=FAR_BASIS, tmode=TM, basis_type=UV, num=6)
    type(ModeInfo), parameter :: MODE_FAR_TE_PQ = ModeInfo(basis=FAR_BASIS, tmode=TE, basis_type=PQ, num=7)
    type(ModeInfo), parameter :: MODE_FAR_TM_PQ = ModeInfo(basis=FAR_BASIS, tmode=TM, basis_type=PQ, num=8)
    type(ModeInfo), parameter :: MODE_BARBER = ModeInfo(basis=BARBER_BASIS, tmode=TETM, basis_type=UV, num=9)
    type(ModeInfo), parameter :: MODE_MISHCH = ModeInfo(basis=MISHCH_BASIS, tmode=TETM, basis_type=UV, num=10)
    type(ModeInfo), parameter :: MODE_FAR_TETM = ModeInfo(basis=FAR_BASIS, tmode = TETM, basis_type=UV, num=11)

    ! type(ModeTransition), parameter :: TRANSITION_SPH_TO_FAR_UV = ModeTransition()

    type :: ModeFactors
        real(knd) :: Qext
        real(knd) :: Qsca
    contains
        procedure :: initialize => initialize_mode_factors
        procedure :: update_and_get_accuracy
        procedure :: qabs => abs_cross_section
    end type ModeFactors

    type :: ModeCalculationResult
        type(ModeFactors) :: factors
        complex(knd), allocatable :: tmatrix(:,:), solution(:)
    contains
        procedure :: initialize => initialize_mode_calculation_result
        final :: delete_mode_calculation_result
    end type ModeCalculationResult

    type :: ScatteringResult
        type(ModeFactors) :: sph_tm, sph_te, far_te, far_tm
    contains
        procedure :: initialize => initialize_scattering_result
        procedure :: update_and_get_accuracy => update_and_get_accuracy_scattering_result
    end type ScatteringResult

    ! type :: ModeQuery
    !     logical :: calculate_factors
    ! contains
    !     procedure :: need => need_calculate_scattering
    ! end type ModeQuery

    ! type :: ScatteringQuery
    !     type(ModeQuery) :: by_mode(SPHEROIDAL_BASIS:MISHCH_BASIS, TE:TM, UV:PQ)
    ! contains
    !     procedure :: need => need_calculate_mode
    ! end type ScatteringQuery

    interface operator(==)
        module procedure mode_info_eq, equal_pairs
    end interface

contains
    character(64) function print_mode_item(this) result(res)
        class(ModeItem), intent(in) :: this

        write(res, *) 'm = ', this%m, 'lnum = ', this%lnum
    end function print_mode_item

    integer function get_matrix_size_by_node(this)
        class(Node), intent(in) :: this

        if (this%info%basis_type == PQ) then
            get_matrix_size_by_node = this%item%lnum
        else
            get_matrix_size_by_node = this%item%lnum * 2
        endif
    end function get_matrix_size_by_node

    type(Node) function construct_node(info, m, lnum, previous, need_calc, to_res) result(res)
        type(ModeInfo), intent(in) :: info
        integer, intent(in) :: m, lnum, previous(:)
        logical, intent(in) :: need_calc, to_res

        res = Node(info, ModeItem(m, lnum), previous, need_calc, to_res)
    end function construct_node

    logical function mode_operator_in(element, list) result(is_in)
        type(ModeInfo), intent(in) :: element, list(1:)
        integer :: i

        is_in = .false.
        do i = 1, size(list)
            if (element == list(i)) then
                is_in = .true.
            endif
        enddo
    end function mode_operator_in

    subroutine delete_node(this)
        type(Node), intent(inout) :: this

        if (allocated(this%previous)) deallocate(this%previous)
    
    end subroutine delete_node

    subroutine initialize_mode_calculation_result(this, ts)
        class(ModeCalculationResult), intent(out) :: this
        integer, intent(in) :: ts

        call this%factors%initialize()

        if (allocated(this%tmatrix) .and. (size(this%tmatrix(:,1)) /= ts)) then
            deallocate(this%tmatrix, this%solution)
        endif

        if (.not. allocated(this%tmatrix)) then
            allocate(this%tmatrix(ts, ts), this%solution(ts))
        endif
    end subroutine initialize_mode_calculation_result

    subroutine delete_mode_calculation_result(this)
        type(ModeCalculationResult), intent(inout) :: this

        if (allocated(this%tmatrix)) then
            deallocate(this%tmatrix) 
        endif

        if (allocated(this%solution)) then
            deallocate(this%solution) 
        endif

    end subroutine delete_mode_calculation_result

    subroutine assert(expression, message)
        logical, intent(in) :: expression
        character(*), intent(in) :: message

        if (.not. expression) then
            write(LOG_FD,*) message
            write(*,*) 'assert failed: ', message, ' see detailes in scattering.log'
            close(LOG_FD)
            call exit(1)
        endif

    end subroutine assert

    function abs_cross_section(this) result(res)
        class(ModeFactors), intent(in) :: this
        real(knd) :: res

        res = this%Qext - this%Qsca

    end function abs_cross_section

    logical function mode_info_eq(this, other)
        class(ModeInfo), intent(in) :: this, other

        mode_info_eq = this%basis == other%basis .and. this%tmode == other%tmode .and. this%basis_type == other%basis_type

    end function mode_info_eq

    logical function equal_pairs(first, second) result(eq)
    type(ModeInfo), intent(in) :: first(2), second(2)

        eq = ((first(1) == second(1)) .and. (first(2) == second(2))) .or. &
            ((first(1) == second(2)) .and. (first(2) == second(1)))
    end function equal_pairs

    function mode_info_to_string(mode) result(res)
        class(ModeInfo), intent(in) :: mode
        character(128) :: res, basis, tmode, basis_type

        if (mode%basis == SPHEROIDAL_BASIS) then
            basis = 'SPH'
        elseif (mode%basis == FAR_BASIS) then
            basis = 'FAR'
        elseif (mode%basis == BARBER_BASIS) then
            basis = 'BARBER'
        elseif (mode%basis == MISHCH_BASIS) then
            basis = 'MISHCH'
        else
            write(*,*) 'logical error: unknown basis'
            basis = 'unknown'
        endif

        if (mode%tmode == TE) then
            tmode = 'TE'
        elseif (mode%tmode == TM) then
            tmode = 'TM'
        elseif (mode%tmode == TETM) then
            tmode = 'TETM'
        else
            write(*,*) 'logical error: unknown tmode'
            tmode = 'unknown'
        endif

        if (mode%basis_type == UV) then
            basis_type = 'UV'
        elseif (mode%basis_type == PQ) then
            basis_type = 'PQ'
        else
            write(*,*) 'logical error: unknown basis_type'
            basis_type = 'unknown'
        endif

        res = trim(basis) // '_'//trim(tmode)//'_'//trim(basis_type)
        ! res = '(' // trim(basis) // ',' // trim(tmode) // ',' trim(basis_type) // ')'
    end function mode_info_to_string

    subroutine initialize_mode_factors(this)
        class(ModeFactors), intent(out) :: this

        this%Qext = 0
        this%Qsca = 0
    
    end subroutine initialize_mode_factors

    subroutine log_mode_factors(name, this)
        type(ModeFactors), intent(in) :: this
        character(*), intent(in) :: name

        write(*,*) name
        write(*,*) 'ext = ', this%Qext
        write(*,*) 'sca = ', this%Qsca
        write(*,*) 'abs = ', this%qabs()

    end subroutine log_mode_factors

    subroutine initialize_scattering_result(this)
        class(ScatteringResult), intent(out) :: this

        call this%sph_te%initialize()
        call this%sph_tm%initialize()
        call this%far_te%initialize()
        call this%far_tm%initialize()
    
    end subroutine initialize_scattering_result

    function update_and_get_accuracy(this, update) result(res)
        class(ModeFactors), intent(inout) :: this
        class(ModeFactors), intent(in) :: update
        real(knd) :: res

        this%Qext = this%Qext + update%Qext
        this%Qsca = this%Qsca + update%Qsca

        if (max(abs(this%Qext), abs(this%Qsca)) < LOWEST_UPDATE) then
            res = LOWEST_UPDATE
        else
            res = max(abs(update%Qext / this%Qext), abs(update%Qsca / this%Qsca))
        endif
    end function update_and_get_accuracy

    function update_and_get_accuracy_scattering_result(this, mode, update) result(res)
        class(ScatteringResult), intent(inout) :: this
        type(ModeInfo), intent(in) :: mode
        class(ModeFactors), intent(in) :: update
        real(knd) :: res

        if (mode == MODE_SPH_TE_PQ .or. mode == MODE_SPH_TE_UV) then
            res = this%sph_te%update_and_get_accuracy(update)
        elseif (mode == MODE_SPH_TM_PQ .or. mode == MODE_SPH_TM_UV) then
            res = this%sph_tm%update_and_get_accuracy(update)
        elseif (mode == MODE_FAR_TE_PQ .or. mode == MODE_FAR_TE_UV) then
            res = this%far_te%update_and_get_accuracy(update)
        elseif (mode == MODE_FAR_TM_PQ .or. mode == MODE_FAR_TM_UV) then
            res = this%far_tm%update_and_get_accuracy(update)
        else
            res = 0
        endif

    end function update_and_get_accuracy_scattering_result

    ! logical function need_calculate_scattering(this)
    !     class(ModeQuery), intent(in) :: this

    !     need_calculate_scattering = this%calculate_factors
    ! end function need_calculate_scattering

    ! logical function need_calculate_mode(this, basis, tmode, potentials) result(res)
    !     class(ScatteringQuery), intent(in) :: this
    !     integer, intent(in) :: basis, tmode, potentials

    !     integer :: i
        
    !     res = this%by_mode(basis, tmode, potentials)%need()

    !     do i = basis + 1, MISHCH_BASIS
    !         res = res .or. this%by_mode(i, tmode, potentials)%need()
    !     enddo
    ! end function need_calculate_mode

    subroutine read_input(filename, f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, model, &
        ntheta, theta0, theta1, nphi, phi0, phi1)
        character(*), intent(in) :: filename
        integer, parameter :: FD_GENERAL = 15
        integer, intent(inout) :: f, nol, matrix_size, spherical_lnum, minm, maxm, ntheta, nphi
        real(knd), allocatable, intent(inout) :: rv(:), xv(:), ab(:)
        real(knd), intent(inout) :: alpha, lambda, theta0, theta1, phi0, phi1
        complex(knd), allocatable, intent(inout) :: ri(:)
        character(32) :: model
    
        character(1024) :: line, name, eq
        integer :: ios, i
        logical :: found_rv, found_xv
    
        open(FD_GENERAL, file=filename)
        read(FD_GENERAL, *) f
        read(FD_GENERAL, *) nol
        allocate(rv(nol), xv(nol), ab(nol), ri(0:nol))
        read(FD_GENERAL, *) xv
        read(FD_GENERAL, *) ab
        read(FD_GENERAL, *) ri
        read(FD_GENERAL, *) lambda
        read(FD_GENERAL, *) alpha
        read(FD_GENERAL, *) matrix_size
        read(FD_GENERAL, *) spherical_lnum
        read(FD_GENERAL, *) minm
        read(FD_GENERAL, *) maxm
        read(FD_GENERAL, *) model
        read(FD_GENERAL, *) ntheta, theta0, theta1
        read(FD_GENERAL, *) nphi, phi0, phi1
    
        close(FD_GENERAL)
    
        if (matrix_size == 0) then
            matrix_size = size_of_matrices(f, xv(1), ab(1), ri(1))
            write(*,*) 'lnum was 0, so it was set to ', matrix_size
        end if
        if (spherical_lnum == 0) then
            spherical_lnum = matrix_size
            write(*,*) 'spherical_lnum was 0, so it was set to ', spherical_lnum
        endif
        alpha = alpha / 180_knd * PI
        theta0 = theta0 / 180_knd * PI
        theta1 = theta1 / 180_knd * PI
        phi0 = phi0 / 180_knd * PI
        phi1 = phi1 / 180_knd * PI
        alpha = max(alpha, 1q-12)
        alpha = min(alpha, PI / 2q0 - 1q-24)
        rv = xv * lambda / (2.0_knd * PI)
    
        write(*, *) 'Read input:'
        write(*, *) 'f = ', f
        write(*, *) 'rv = ', rv
        write(*, *) 'xv = ', xv
        write(*, *) 'ab = ', ab
        write(*, *) 'alpha = ', alpha, 'radian = ', alpha / PI * 180_knd, 'degrees'
        write(*, *) 'lambda = ', lambda
        write(*, *) 'ri = ', ri
        write(*, *) 'lnum = ', matrix_size
        write(*, *) 'spherical_lnum = ', spherical_lnum
        write(*, *) 'm = ', minm, ':', maxm
        write(*,*) 'model = ', model
    end subroutine read_input
end module utils