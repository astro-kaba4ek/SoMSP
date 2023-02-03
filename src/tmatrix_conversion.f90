module tmatrix_conversion
    use regime
    use constants
    use spheroidal
    use legendre_functions
    use contexts
    use utils
    implicit none

    private
    public :: convert_tmatrix

    type, abstract :: ModeTransition
    contains
        procedure(set_new_tmatrix_from_old), deferred, nopass :: convert_tmatrix
    end type ModeTransition

    abstract interface
        subroutine set_new_tmatrix_from_old(&
                computation_context, &
                old_item, old_tmatrix, &
                new_item, new_tmatrix)
            import :: knd
            import :: ComputationContext, ModeItem

            type(ComputationContext), intent(inout) :: computation_context
            type(ModeItem), intent(in) :: old_item, new_item
            complex(knd), intent(in) :: old_tmatrix(:,:)
            complex(knd), intent(out) :: new_tmatrix(:,:)

        end subroutine set_new_tmatrix_from_old
    end interface

    type, extends(ModeTransition) :: SphToFarUVTransition
    contains
        procedure, nopass :: convert_tmatrix => set_far_spherical_from_spheroidal_uv
    end type SphToFarUVTransition

    type, extends(ModeTransition) :: FarToSphUVTransition
    contains
        procedure, nopass :: convert_tmatrix => set_spheroidal_from_far_spherical_uv
    end type FarToSphUVTransition


    type, extends(ModeTransition) :: SphToFarPQTransition
    contains
        procedure, nopass :: convert_tmatrix => set_far_from_sph_pq
    end type SphToFarPQTransition

    type, extends(ModeTransition) :: FarToSphPQTransition
    contains
        procedure, nopass :: convert_tmatrix => set_spheroidal_from_far_spherical_pq
    end type FarToSphPQTransition


    type, extends(ModeTransition) :: FarUVTEToBarberTransition
    contains
        procedure, nopass :: convert_tmatrix => set_barber_from_far_uv_te
    end type FarUVTEToBarberTransition

    type, extends(ModeTransition) :: BarberToFarUVTETransition
    contains
        procedure, nopass :: convert_tmatrix => set_far_uv_te_from_barber
    end type BarberToFarUVTETransition


    type, extends(ModeTransition) :: FarUVTMToBarberTransition
    contains
        procedure, nopass :: convert_tmatrix => set_barber_from_far_uv_tm
    end type FarUVTMToBarberTransition

    type, extends(ModeTransition) :: BarberToFarUVTMTransition
    contains
        procedure, nopass :: convert_tmatrix => set_far_uv_tm_from_barber
    end type BarberToFarUVTMTransition


    type, extends(ModeTransition) :: FarPQTETMToBarberTransition
    contains
        procedure, nopass :: convert_tmatrix => set_barber_from_far_pq_tetm
    end type FarPQTETMToBarberTransition


    type, extends(ModeTransition) :: BarberToMishchTransition
    contains
        procedure, nopass :: convert_tmatrix => set_mishch_from_barber
    end type BarberToMishchTransition

    type, extends(ModeTransition) :: MishchToBarberTransition
    contains
        procedure, nopass :: convert_tmatrix => set_barber_from_mishch
    end type MishchToBarberTransition

contains

    subroutine convert_tmatrix(computation_context, old_node, new_node, old_tmatrix, new_tmatrix)
        type(Node), intent(in) :: old_node, new_node
        type(ComputationContext), intent(inout) :: computation_context
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(ModeInfo) :: old_mode, new_mode

        integer :: lnum, spherical_lnum

        class(ModeTransition), allocatable :: transition

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} convert tmatrix'
        transition = build_transition(old_node%info, new_node%info)

        if (LOG_INFO) write(LOG_FD, *) '{INFO} convert tmatrix from '//trim(old_node%info%to_string())//' {'&
        //trim(old_node%item%to_string())&
        //'} to '//trim(new_node%info%to_string())//' {'//trim(new_node%item%to_string())//'}'

        call transition%convert_tmatrix(computation_context, old_node%item, old_tmatrix, new_node%item, new_tmatrix)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} convert tmatrix'

        ! deallocate(transition)

    end subroutine convert_tmatrix

    function build_transition(old_mode, new_mode) result(tr)
        type(ModeInfo), intent(in) :: old_mode, new_mode
        class(ModeTransition), allocatable :: tr

        if (allocated(tr)) deallocate(tr)
        ! sph <-> far : uv
        if     ((old_mode == MODE_SPH_TE_UV .and. new_mode == MODE_FAR_TE_UV) .or. &
                (old_mode == MODE_SPH_TM_UV .and. new_mode == MODE_FAR_TM_UV)) then
            allocate(SphToFarUVTransition::tr)
        elseif ((old_mode == MODE_FAR_TE_UV .and. new_mode == MODE_SPH_TE_UV) .or. &
                (old_mode == MODE_FAR_TM_UV .and. new_mode == MODE_SPH_TM_UV)) then
            allocate(FarToSphUVTransition::tr)
        ! sph <-> far : pq
        elseif ((old_mode == MODE_SPH_TE_PQ .and. new_mode == MODE_FAR_TE_PQ) .or. &
                (old_mode == MODE_SPH_TM_PQ .and. new_mode == MODE_FAR_TM_PQ)) then
            allocate(SphToFarPQTransition::tr)
        elseif ((old_mode == MODE_FAR_TE_PQ .and. new_mode == MODE_SPH_TE_PQ) .or. &
                (old_mode == MODE_FAR_TM_PQ .and. new_mode == MODE_SPH_TM_PQ)) then
            allocate(FarToSphPQTransition::tr)           
        ! far uv te <-> barber
        elseif (old_mode == MODE_FAR_TE_UV .and. new_mode == MODE_BARBER) then
            allocate(FarUVTEToBarberTransition::tr)
        elseif (old_mode == MODE_BARBER .and. new_mode == MODE_FAR_TE_UV) then
            allocate(BarberToFarUVTETransition::tr)
        ! far uv tm <-> barber
        elseif (old_mode == MODE_FAR_TM_UV .and. new_mode == MODE_BARBER) then
            allocate(FarUVTMToBarberTransition::tr)
        elseif (old_mode == MODE_BARBER .and. new_mode == MODE_FAR_TM_UV) then
            allocate(BarberToFarUVTMTransition::tr)
        ! far pq tetm -> barber
        elseif (old_mode == MODE_FAR_TETM .and. new_mode == MODE_BARBER) then
            allocate(FarPQTETMToBarberTransition::tr)
        ! barber <-> mishch
        elseif (old_mode == MODE_BARBER .and. new_mode == MODE_MISHCH) then
            allocate(BarberToMishchTransition::tr)
        elseif (old_mode == MODE_MISHCH .and. new_mode == MODE_BARBER) then
            allocate(MishchToBarberTransition::tr)
        else
            call assert(.false., 'No know conversion between modes '//trim(old_mode%to_string())// &
            ' and '//trim(new_mode%to_string()))
        endif

    end function build_transition

    subroutine multiply_with_transpose(spheroidal_lnum, spheroidal_tmatrix, connecting_matrix, spherical_lnum, spherical_tmatrix)
        integer, intent(in) :: spheroidal_lnum, spherical_lnum
        complex(knd), intent(in) :: spheroidal_tmatrix(spheroidal_lnum, spheroidal_lnum), &
                connecting_matrix(spheroidal_lnum, spherical_lnum)
        complex(knd), intent(out) :: spherical_tmatrix(spherical_lnum, spherical_lnum)

        spherical_tmatrix = matmul(transpose(connecting_matrix), matmul(spheroidal_tmatrix, connecting_matrix))

    end subroutine multiply_with_transpose

    subroutine set_far_from_sph_pq(computation_context, &
        old_item, old_tmatrix, &
        new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(FarContext), pointer :: context

        call check_valid_pq_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'sph_pq -> far_pq')

        context => computation_context%get_far_context(old_item%m, old_item%lnum, new_item%lnum)

        call multiply_with_transpose( &
            old_item%lnum, &
            old_tmatrix, &
            context%connecting_matrix, &
            new_item%lnum, &
            new_tmatrix &
        )

    end subroutine set_far_from_sph_pq

    subroutine set_spheroidal_from_far_spherical_pq(computation_context, &
        old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: new_item, old_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(FarContext), pointer :: context

        call check_valid_pq_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'far_pq -> sph_pq')

        context => computation_context%get_far_context(old_item%m, new_item%lnum, old_item%lnum)

        call multiply_with_transpose( &
            old_item%lnum, &
            old_tmatrix, &
            transpose(context%connecting_matrix),&
            new_item%lnum, &
            new_tmatrix &
        )

    end subroutine set_spheroidal_from_far_spherical_pq

    subroutine set_far_spherical_from_spheroidal_uv(computation_context, &
        old_item, old_tmatrix, &
        new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(FarContext), pointer :: context

        integer :: i, j

        call check_valid_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'sph_uv -> far_uv')

        context => computation_context%get_far_context(old_item%m, old_item%lnum, new_item%lnum)

        do i = 0, 1
            do j = 0, 1
                call multiply_with_transpose( &
                    old_item%lnum, &
                    old_tmatrix((i * old_item%lnum + 1):((i + 1) * old_item%lnum), &
                                       (j * old_item%lnum + 1):((j + 1) * old_item%lnum)), &
                    context%connecting_matrix,&
                    new_item%lnum, &
                    new_tmatrix((i * new_item%lnum + 1):((i + 1) * new_item%lnum), &
                                      (j * new_item%lnum + 1):((j + 1) * new_item%lnum)) &
                )
            end do
        end do

    end subroutine set_far_spherical_from_spheroidal_uv

    subroutine set_spheroidal_from_far_spherical_uv(computation_context, &
        old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: new_item, old_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        type(FarContext), pointer :: context

        integer :: i, j

        call check_valid_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'far_uv -> sph_uv')

        context => computation_context%get_far_context(old_item%m, new_item%lnum, old_item%lnum)

        do i = 0, 1
            do j = 0, 1
                call multiply_with_transpose( &
                    old_item%lnum, &
                    old_tmatrix((i * old_item%lnum + 1):((i + 1) * old_item%lnum), &
                                    (j * old_item%lnum + 1):((j + 1) * old_item%lnum)), &
                    transpose(context%connecting_matrix),&
                    new_item%lnum, &
                    new_tmatrix((i * new_item%lnum + 1):((i + 1) * new_item%lnum), &
                                       (j * new_item%lnum + 1):((j + 1) * new_item%lnum)) &
                )
            end do
        end do

    end subroutine set_spheroidal_from_far_spherical_uv

    subroutine set_barber_from_far_uv_te(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :)
        complex(knd), intent(out) :: new_tmatrix(:, :)

        complex(knd) :: orth(old_item%lnum, old_item%lnum)
        integer :: m, lnum, i, j, i1, j1
        real(knd) :: k
        type(BarberContext), pointer :: context

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'far_uv_te -> barber')

        m = old_item%m
        lnum = old_item%lnum
        k = computation_context%scattering_context%calculation_point%k

        context => computation_context%get_barber_context(m, lnum)

        ! for T_VV^{11}
        orth = old_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) + k * get_double_mult_left(context%f, lnum,&
        old_tmatrix(1:lnum, (lnum + 1):(2 * lnum)))
        call multiply_by_diag_right(orth, lnum, context%xd)
        call multiply_by_diag_left(orth, lnum, context%xdrev)
        new_tmatrix(1:lnum, 1:lnum) = orth
        
        ! for T_VV^{12}
        orth = old_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) + k * get_double_mult_left(context%f, lnum, &
        old_tmatrix(1:lnum, (lnum + 1):(2 * lnum)))
        orth = get_double_mult_right(orth, lnum, context%f)
        orth = get_double_mult_left(context%f, lnum, old_tmatrix(1:lnum,1:lnum)) + &
        1.0_knd / k * old_tmatrix((lnum + 1):(2 * lnum), 1:lnum) - orth
        call multiply_by_diag_right(orth, lnum, context%xd * context%grev)
        call multiply_by_diag_left(orth, lnum, context%xdrev)
        new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) = orth

        ! for T_{VV}^{21}
        orth = old_tmatrix(1:lnum, (lnum + 1):(2 *lnum))
        call multiply_by_diag_right(orth, lnum, context%xd)
        call multiply_by_diag_left(orth, lnum, k * context%g * context%xdrev)
        new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) = orth

        ! for T_{VV}^{22}
        orth = old_tmatrix(1:lnum, 1:lnum) - k * get_double_mult_right(old_tmatrix(1:lnum, (lnum + 1):(2 * lnum)), lnum, context%f)
        call multiply_by_diag_right(orth, lnum, context%xd * context%grev)
        call multiply_by_diag_left(orth, lnum, context%xdrev * context%g)
        new_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) = orth

    end subroutine set_barber_from_far_uv_te
    
    function ff(x)
        real :: x, ff
    end function ff

    function get_double_mult_right(a, n, b) result(c)
        complex(knd) :: a(n,n), b(n,n), c(n,n)
        integer :: n, i, j
        c = 0
        do j = 1, n
            do i = 1, n
                if (j > 1) then
                    c(i,j) = c(i,j) + a(i, j - 1) * b(j - 1, j)
                endif
                if (j < n) then
                    c(i,j) = c(i,j) + a(i, j + 1) * b(j + 1, j)
                endif
            enddo
        enddo
        
    end function get_double_mult_right

    function get_double_mult_left(a, n, b) result(c)
        complex(knd) :: a(n,n), b(n,n), c(n,n)
        integer :: n, i, j
        c = 0
        do j = 1, n
            do i = 1, n
                if (i > 1) then
                    c(i,j) = c(i,j) + a(i, i - 1) * b(i - 1, j)
                endif
                if (i < n) then
                    c(i,j) = c(i,j) + a(i, i + 1) * b(i + 1, j)
                endif
            enddo
        enddo
        
    end function get_double_mult_left

    subroutine set_barber_from_far_uv_tm(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        integer :: lnum, i, j, i1, j1

        complex(knd) :: intermediate(old_item%lnum, old_item%lnum)

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'far_uv_tm -> barber')

        call set_barber_from_far_uv_te(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        
        lnum = old_item%lnum

        intermediate = new_tmatrix(1:lnum, (lnum + 1):(2 * lnum))
        new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) = -new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) 
        new_tmatrix((lnum + 1):(2 * lnum), 1:lnum)  = -intermediate

        intermediate = new_tmatrix(1:lnum, 1:lnum) 
        new_tmatrix(1:lnum, 1:lnum) = new_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) 
        new_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum))  = intermediate

    end subroutine set_barber_from_far_uv_tm

    subroutine set_far_uv_te_from_barber(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :)
        complex(knd), intent(out) :: new_tmatrix(:, :)

        type(BarberContext), pointer :: context

        complex(knd) :: p(2 * old_item%lnum, 2 * old_item%lnum), tmp(old_item%lnum, old_item%lnum)
        integer :: m, lnum, i, j, i1, j1
        real(knd) :: k

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'barber -> far_uv_te')

        m = old_item%m
        lnum = old_item%lnum
        k = computation_context%scattering_context%calculation_point%k

        context => computation_context%get_barber_context(m, lnum)
 
        p = old_tmatrix
        do i = 0, 1
            do j = 0, 1
                call multiply_by_diag_right(p((i * lnum + 1):(i + 1) * lnum, (j * lnum + 1):(j + 1) * lnum), lnum, &
                    context%xdrev)
                call multiply_by_diag_left(p((i * lnum + 1):(i + 1) * lnum, (j * lnum + 1):(j + 1) * lnum), lnum, &
                    context%xd)
            enddo
        enddo
        
        ! for T_UV^{11}
        tmp = p((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum))
        call multiply_by_diag_right(tmp, lnum, context%g)
        tmp = tmp + get_double_mult_right(p((lnum + 1):(2 * lnum), 1:lnum), lnum, context%f)
        call multiply_by_diag_left(tmp, lnum, context%grev)
        new_tmatrix(1:lnum, 1:lnum) = tmp

        ! for T_UV^{12}

        new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) = p((lnum + 1):(2 * lnum), 1:lnum)
        call multiply_by_diag_left(new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)), lnum, context%grev / k)

        ! for T_{UV}^{21}
        tmp = p(1:lnum, (lnum + 1):(2 * lnum))
        call multiply_by_diag_right(tmp, lnum, context%g)
        new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) = tmp + get_double_mult_right(p(1:lnum, 1:lnum), lnum, context%f) -&
            get_double_mult_left(context%f, lnum, new_tmatrix(1:lnum, 1:lnum))

        ! for T_{UV}^{22}
        tmp = p((lnum + 1):(2 * lnum), 1:lnum)
        call multiply_by_diag_left(tmp, lnum, context%grev)
        new_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) = p(1:lnum, 1:lnum) - &
            get_double_mult_left(context%f, lnum, tmp)

    end subroutine set_far_uv_te_from_barber

    subroutine set_far_uv_tm_from_barber(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :)
        complex(knd), intent(out) :: new_tmatrix(:, :)

        complex(knd) :: intermediate(2 * old_item%lnum, 2 * old_item%lnum)
        integer :: lnum

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'barber -> far_uv_tm')

        lnum = old_item%lnum

        intermediate(1:lnum, 1:lnum) = old_tmatrix((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum))
        intermediate((lnum + 1):(2 * lnum), (lnum + 1):(2 * lnum)) = old_tmatrix(1:lnum, 1:lnum)
        intermediate(1:lnum, (lnum + 1):(2 * lnum)) = -old_tmatrix((lnum + 1):(2 * lnum), 1:lnum)
        intermediate((lnum + 1):(2 * lnum), 1:lnum) = -old_tmatrix(1:lnum, (lnum + 1):(2 * lnum))

        call set_far_uv_te_from_barber(computation_context, old_item, intermediate, new_item, new_tmatrix)

    end subroutine set_far_uv_tm_from_barber

    subroutine set_barber_from_far_pq_tetm(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:,:)
        complex(knd), intent(out) :: new_tmatrix(:,:)

        integer :: lnum, i, j, i1, j1

        complex(knd) :: h(old_item%lnum)

        lnum = old_item%lnum
        call assert(old_item%lnum == new_item%lnum, 'different lnum in far_pq_tetm -> barber')

        call assert(old_item%m == 1, 'wrong far m in far_pq_tetm -> barber')
        call assert(new_item%m == 0, 'wrong barber m in far_pq_tetm -> barber')

        call assert(all(shape(old_tmatrix) == [2 * lnum, 2 * lnum]), &
        'wrong old tmatrix shape in far_pq_tetm -> barber')
        call assert(all(shape(new_tmatrix) == [2 * lnum, 2 * lnum]), &
        'wrong new tmatrix shape in far_pq_tetm -> barber')

        new_tmatrix = -old_tmatrix

        h = [(sqrt(i * (i + 1q0) / (2q0 * i + 1q0)), i = 1, lnum)]

        write(LOG_FD,*) 'h = ', h
        write(LOG_FD,*) 'h-1 = ', cmplx(1q0, 0q0) / h

        do i = 0, 1
            call multiply_by_diag_right(new_tmatrix(i * lnum + 1:(i + 1)*lnum, i * lnum + 1:(i + 1)*lnum), &
            lnum, cmplx(1q0, 0q0) / h)
            call multiply_by_diag_left(new_tmatrix(i * lnum + 1:(i + 1)*lnum, i * lnum + 1:(i + 1)*lnum), &
            lnum, cmplx(1q0, 0q0) * (h))
        enddo

    end subroutine set_barber_from_far_pq_tetm
    

    subroutine set_mishch_from_barber(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :)
        complex(knd), intent(out) :: new_tmatrix(:, :)

        type(BarberContext), pointer :: context

        integer :: m, lnum, i, j, i1, j1

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'barber -> mishch')

        m = old_item%m
        lnum = old_item%lnum
        context => computation_context%get_barber_context(m, lnum)

        new_tmatrix = old_tmatrix

        do i1 = 0, 1
            do j1 = 0, 1
                call multiply_by_diag_right(new_tmatrix((i1 * lnum + 1):(i1 + 1) * lnum, (j1 * lnum + 1):(j1 + 1) * lnum), &
                    lnum, context%mishch_mult)
                call multiply_by_diag_left(new_tmatrix((i1 * lnum + 1):(i1 + 1) * lnum, (j1 * lnum + 1):(j1 + 1) * lnum), &
                    lnum, 1.0_knd / context%mishch_mult)
            end do
        end do
        new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) = -new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) * ideg(1)
        new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) = new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) * ideg(1)

    end subroutine set_mishch_from_barber

    subroutine set_barber_from_mishch(computation_context, old_item, old_tmatrix, new_item, new_tmatrix)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :)
        complex(knd), intent(out) :: new_tmatrix(:, :)

        type(BarberContext), pointer :: context

        integer :: m, lnum, i, j, i1, j1

        call check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, 'mishch -> barber')

        m = new_item%m
        lnum = new_item%lnum

        context => computation_context%get_barber_context(m, lnum)

        new_tmatrix = old_tmatrix

        new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) = new_tmatrix(1:lnum, (lnum + 1):(2 * lnum)) * ideg(1)
        new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) = -new_tmatrix((lnum + 1):(2 * lnum), 1:lnum) * ideg(1)

        do i1 = 0, 1
            do j1 = 0, 1
                call multiply_by_diag_right(new_tmatrix((i1 * lnum + 1):(i1 + 1) * lnum, (j1 * lnum + 1):(j1 + 1) * lnum), &
                    lnum, 1.0_knd / context%mishch_mult)
                call multiply_by_diag_left(new_tmatrix((i1 * lnum + 1):(i1 + 1) * lnum, (j1 * lnum + 1):(j1 + 1) * lnum), &
                    lnum, context%mishch_mult)
            end do
        end do

    end subroutine set_barber_from_mishch

    subroutine check_valid_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, transition_name)
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :), new_tmatrix(:, :)
        character(*), intent(in) :: transition_name

        call assert(old_item%m == new_item%m, 'different m in '//trim(transition_name))

        call assert(all(shape(old_tmatrix) == [2 * old_item%lnum, 2 * old_item%lnum]), &
        'wrong old tmatrix shape in '//trim(transition_name))
        call assert(all(shape(new_tmatrix) == [2 * new_item%lnum, 2 * new_item%lnum]), &
        'wrong new tmatrix shape in '//trim(transition_name))

    end subroutine check_valid_uv_transition

    subroutine check_valid_equal_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, transition_name)
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :), new_tmatrix(:, :)
        character(*), intent(in) :: transition_name

        call assert(old_item%lnum == new_item%lnum, 'different lnum in '//trim(transition_name))

        call check_valid_uv_transition(old_item, old_tmatrix, new_item, new_tmatrix, transition_name)

    end subroutine check_valid_equal_uv_transition

    subroutine check_valid_pq_transition(old_item, old_tmatrix, new_item, new_tmatrix, transition_name)
        type(ModeItem), intent(in) :: old_item, new_item
        complex(knd), intent(in) :: old_tmatrix(:, :), new_tmatrix(:, :)
        character(*), intent(in) :: transition_name

        call assert(old_item%m == 1, 'wrong old m in '//trim(transition_name))
        call assert(new_item%m == 1, 'wrong new m in '//trim(transition_name))

        call assert(all(shape(old_tmatrix) == [old_item%lnum, old_item%lnum]), 'wrong old tmatrix shape in '//trim(transition_name))
        call assert(all(shape(new_tmatrix) == [new_item%lnum, new_item%lnum]), 'wrong new tmatrix shape in '//trim(transition_name))

    end subroutine check_valid_pq_transition

end module tmatrix_conversion