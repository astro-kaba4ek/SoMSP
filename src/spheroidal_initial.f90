! Created by drakosha on 15.07.2021.

module spheroidal_initial
    use regime
    use spheroidal
    use wavelength_point
    use regime
    use angle
    use constants
    use contexts
    implicit none

contains
    ! Initializers
    subroutine set_initial_pq_te(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        integer :: i, lnum
        type(SpheroidalContext), pointer :: context

        lnum = mode_item%lnum

        if (allocated(initial) .and. size(initial) /= lnum) then
            deallocate(initial)
        endif

        if(.not. allocated(initial)) then
            allocate(initial(lnum))
        endif

        context => computation_context%get_spheroidal_context(1, mode_item%m, lnum, lnum)

        do i = 1, lnum
            initial(i) = -2q0 * IDEG(mod(i, 4)) * context%layers(0,1)%s1(i, 1)
        enddo

        ! initial = initial

        !call log_array('outside_layer%s1', matrix_size, outside_layer%s1(:, argnum), 10, FD_INFO)
        !write(*,*) 'lnum = ', outside_layer%lnum
        !call log_array('symmetric_initial_te', matrix_size, initial, 10, FD_INFO)
    end subroutine set_initial_pq_te

    subroutine set_initial_uv_te(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        integer :: i, j, lnum
        real(knd) :: coefficient, alpha, k
        type(SpheroidalContext), pointer :: context

        lnum = mode_item%lnum

        if (allocated(initial) .and. size(initial) /= 2 * lnum) then
            deallocate(initial)
        endif

        if(.not. allocated(initial)) then
            allocate(initial(2*lnum))
        endif

        initial = 0
        k = computation_context%scattering_context%calculation_point%k

        context => computation_context%get_spheroidal_context(1, mode_item%m, lnum, lnum)

        if (abs(computation_context%scattering_context%directions%alpha%value) > 1q-6) then
            coefficient = context%layers(0,1)%spheroidal_type * 4.0q0 / (&
            k * computation_context%scattering_context%directions%alpha%angle_sin)
            do i = 1, lnum
                initial(i) = coefficient * IDEG(mod(i + context%layers(0,1)%m + 2, 4)) * context%layers(0,1)%s1(i, 1)
            enddo
        elseif (context%layers(0,1)%m == 1) then
            coefficient = context%layers(0,1)%spheroidal_type * 4.0q0 / k
            do i = 1, lnum
                do j = 1 - mod(i, 2), context%layers(0,1)%maxd, 2
                    initial(i) = initial(i) + (j + 1) * (j + 2) / 2.0_knd * context%layers(0,1)%legendre(j, i)
                end do
                initial(i) = coefficient * IDEG(mod(i + context%layers(0,1)%m + 2, 4)) * initial(i)
            end do
        endif

        initial = -initial * context%layers(0,1)%spheroidal_type

    end subroutine set_initial_uv_te

    subroutine set_initial_pq_tm(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        integer :: i
        real(knd) :: eps_to_mu
        type(WavelengthPoint), pointer :: wavelength

        wavelength => computation_context%scattering_context%calculation_point
        eps_to_mu = sqrt(wavelength%eps(0) / wavelength%mu(0))

        call set_initial_pq_te(computation_context, mode_item, initial)
        initial = initial * eps_to_mu

    end subroutine set_initial_pq_tm

    subroutine set_initial_uv_tm(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        integer :: i
        real(knd) :: eps_to_mu
        type(WavelengthPoint), pointer :: wavelength

        wavelength => computation_context%scattering_context%calculation_point

        eps_to_mu = sqrt(wavelength%eps(0) / wavelength%mu(0))

        call set_initial_uv_te(computation_context, mode_item, initial)
        initial = initial * eps_to_mu

    end subroutine set_initial_uv_tm

    subroutine set_initial_far_uv_te(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        type(LegendreCalculation) :: legendre
        integer :: m, lnum, i
        type(WavelengthPoint), pointer :: calculation_point
        type(AngleType), pointer :: alpha

        m = mode_item%m
        lnum = mode_item%lnum
        calculation_point => computation_context%scattering_context%calculation_point
        alpha => computation_context%scattering_context%directions%alpha
        call legendre%set(m, lnum + m - 1, alpha%angle_cos)
        call legendre%calculate()


        if (allocated(initial) .and. size(initial) /= 2 * lnum) then
            deallocate(initial)
        endif

        if(.not. allocated(initial)) then
            allocate(initial(2 * lnum))
        endif
        initial = 0

        do i = 1, lnum
            initial(i) = IDEG(mod(i + legendre%m + 2, 4)) * legendre%pr(1, i) * sqrt(legendre%coef(i))
        enddo

        initial = initial * (4q0 / (calculation_point%k * alpha%angle_sin))

    end subroutine set_initial_far_uv_te

    subroutine set_initial_far_pq_tm(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        call set_initial_far_pq_te(computation_context, mode_item, initial)
        initial = initial

    end subroutine set_initial_far_pq_tm

    subroutine set_initial_far_pq_te(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        type(LegendreCalculation) :: legendre
        integer :: m, lnum, i
        type(AngleType), pointer :: alpha

        m = mode_item%m
        lnum = mode_item%lnum
        alpha => computation_context%scattering_context%directions%alpha
        call legendre%set(m, lnum + m - 1, alpha%angle_cos)
        call legendre%calculate()

        if (allocated(initial) .and. size(initial) /= lnum) then
            deallocate(initial)
        endif

        if(.not. allocated(initial)) then
            allocate(initial(lnum))
        endif

        initial = 0

        do i = 1, lnum
            initial(i) = -2q0 * IDEG(mod(i + m + 3, 4)) * legendre%pr(1, i) * sqrt(legendre%coef(i))
        end do

    end subroutine set_initial_far_pq_te

    subroutine set_initial_far_uv_tm(computation_context, mode_item, initial)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), allocatable, intent(out) :: initial(:)

        call set_initial_far_uv_te(computation_context, mode_item, initial)
        initial = initial

    end subroutine set_initial_far_uv_tm

end module spheroidal_initial