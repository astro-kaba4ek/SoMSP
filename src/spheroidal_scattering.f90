! Created by drakosha on 15.07.2021.

module spheroidal_scattering
    use regime
    use spheroidal_tmatrix
    use spheroidal_initial
    use spheroidal
    use spheroidal_scatterer
    use wavelength_point
    use angle
    use legendre_functions
    implicit none
    integer, private, parameter :: general_log_fd = 111
contains

    ! Extinction factor
    function get_extinction_factor_sph_pq(computation_context, mode_item, solution) result(ext)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        real(knd) :: ext

        integer :: i
        type(SpheroidalContext), pointer :: context

        context => computation_context%get_spheroidal_context(1, mode_item%m, mode_item%lnum, mode_item%lnum)
        ext = 0

        do i = 1, mode_item%lnum
            ext = ext + real(solution(i) * IDEG(mod(3 * i, 4)) * context%layers(0,1)%s1(i, 1), knd)
        enddo

        ext = -ext * 4.0q0 * computation_context%scattering_context%scatterer%common_factor(1)

    end function get_extinction_factor_sph_pq

    function get_extinction_factor_sph_uv(computation_context, mode_item, solution) result(ext)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        real(knd) :: ext, k1
        integer :: i, m, lnum
        type(SpheroidalContext), pointer :: context

        ext = 0
        k1 = computation_context%scattering_context%calculation_point%k
        m = mode_item%m
        lnum = mode_item%lnum

        context => computation_context%get_spheroidal_context(1, m, lnum, lnum)

        do i = 1, lnum
            ext = ext + real(NEGIDEG(mod(i + m + 2, 4)) * &
                (k1 * solution(i) * context%layers(0,1)%s1(i, 1) - &
                    solution(i + lnum) * NEGIDEG(1) * context%layers(0,1)%s1d(i, 1)), knd)
        enddo

        if (m == 0) then
            ext = ext * 0.5_knd
        end if

        ext = ext * 4.0_knd * &
            computation_context%scattering_context%scatterer%common_factor(1) * &
            computation_context%scattering_context%directions%alpha%angle_sin

    end function get_extinction_factor_sph_uv

    function get_scattering_factor_sph_pq(computation_context, mode_item, solution)result(sca)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        real(knd) :: sca

        sca = sum(abs(solution)**2) * 2.0_knd * computation_context%scattering_context%scatterer%common_factor(1)

    end function get_scattering_factor_sph_pq

    function get_scattering_factor_sph_uv(computation_context, mode_item, solution) result(sca)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        integer :: i, j, m, md, lnum
        complex(knd), allocatable, dimension(:,:) :: Omega, Kappa, Tau
        real(knd) :: sca, k1
        real(knd), allocatable :: mult_coef(:)
        type(SpheroidalContext), pointer :: context

        k1 = computation_context%scattering_context%calculation_point%k
        sca = 0

        m = mode_item%m
        lnum = mode_item%lnum
        context => computation_context%get_spheroidal_context(1, m, lnum, lnum)

        allocate(Omega(lnum, lnum), Kappa(lnum, lnum), Tau(lnum, lnum))

        md = context%layers(0,1)%maxd

        call fill_common_multiplier(m, md, mult_coef)

        ! write(*,*) 'start omega calculation'
        call triple_dep_integral(Omega, m, context%layers(0,1)%legendre, context%layers(0,1)%legendre, &
        mult_coef, omega_c_lower, &
        omega_c_middle, &
        omega_c_upper)
        ! write(*,*) 'end omega calculation'
        ! call calculate_omega(outside_layer, outside_layer, Omega, matrix_size)
        !call calculate_kappa(outside_layer, outside_layer, Kappa, matrix_size)
        call double_dep_integral(Kappa, m, context%layers(0,1)%legendre, context%layers(0,1)%legendre, &
        mult_coef, kappa_c_lower, kappa_c_upper)
        ! call calculate_tau(outside_layer, outside_layer, Tau, matrix_size)
        call single_dep_integral(Tau, m, context%layers(0,1)%legendre, context%layers(0,1)%legendre, &
        mult_coef, tau_c_middle)
        deallocate(mult_coef)

        do i = 1, lnum
            do j = 1, lnum
                sca = sca + real(IDEG(mod(j + i * 3, 4)) * (&
                        k1**2 * solution(i) * conjg(solution(j)) * Omega(i, j) + &
                                IDEG(1) * k1 * (&
                                        solution(i + lnum) * conjg(solution(j)) * Kappa(j, i) - &
                                                conjg(solution(j + lnum)) * solution(i) * Kappa(i, j)) + &
                                solution(i + lnum) * conjg(solution(j + lnum)) * &
                                        Tau(i, j)), knd)
            end do
        end do
        if (m == 0) then
            sca = sca * 0.5_knd
        end if

        sca = sca * computation_context%scattering_context%scatterer%common_factor(1)

        deallocate(Omega, Kappa, Tau)

    end function get_scattering_factor_sph_uv

    real(knd) function get_extinction_factor_far_pq(computation_context, mode_item, solution) result(ext)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        type(LegendreCalculation) :: legendre
        integer :: m, lnum, i
        type(WavelengthPoint), pointer :: calculation_point
        type(AngleType), pointer :: alpha
        type(SpheroidalShape), pointer :: shape

        m = mode_item%m
        lnum = mode_item%lnum
        calculation_point => computation_context%scattering_context%calculation_point
        alpha => computation_context%scattering_context%directions%alpha
        call legendre%set(m, lnum + m - 1, alpha%angle_cos)
        call legendre%calculate()

        ext = 0q0

        do i = 1, lnum
            ! write(*,*) NEGIDEG(mod(i + legendre%m + 3, 4)), solution(i) , legendre%pr(1, i), sqrt(legendre%coef(i))
            ext = ext + real(NEGIDEG(mod(i + m + 3, 4)) * solution(i) * legendre%pr(1, i) * sqrt(legendre%coef(i)), knd)
        enddo

        shape => computation_context%scattering_context%scatterer%shape
        ext = convertCtoQ(-ext * 4q0 * PI / calculation_point%k ** 2, shape%spheroidal_type, shape%rv, shape%ab, shape%alpha%value)
        ! write(*,*) 'ext2 = ', ext
    end function get_extinction_factor_far_pq

    real(knd) function get_extinction_factor_far_uv(computation_context, mode_item, solution) result(ext)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        type(LegendreCalculation) :: legendre
        integer :: m, lnum, i
        type(WavelengthPoint), pointer :: calculation_point
        type(AngleType), pointer :: alpha
        type(SpheroidalShape), pointer :: shape

        m = mode_item%m
        lnum = mode_item%lnum
        calculation_point => computation_context%scattering_context%calculation_point
        alpha => computation_context%scattering_context%directions%alpha
        call legendre%set(m, lnum + m - 1, alpha%angle_cos)
        call legendre%calculate()

        ext = 0q0

        do i = 1, lnum
            ext = ext + real(NEGIDEG(mod(i + m + 2, 4)) * (calculation_point%k * solution(i) * legendre%pr(1, i) + &
                    solution(i + lnum) * IDEG(1) * legendre%pdr(1, i)), knd) * sqrt(legendre%coef(i))! / (this%legendre%coef(i))
        enddo

        shape => computation_context%scattering_context%scatterer%shape
        ext = convertCtoQ(-ext * 4q0 * PI / calculation_point%k ** 2 * alpha%angle_sin, &
        shape%spheroidal_type, shape%rv, shape%ab, shape%alpha%value)
        if (m == 0) then
            ext = ext / 2
        end if
    end function get_extinction_factor_far_uv

    real(knd) function get_scattering_factor_far_pq(computation_context, mode_item, solution) result(sca)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        integer :: lnum, i
        type(WavelengthPoint), pointer :: calculation_point
        type(SpheroidalShape), pointer :: shape

        lnum = mode_item%lnum
        calculation_point => computation_context%scattering_context%calculation_point

        sca = 0

        do i = 1, lnum
            sca = sca + abs(solution(i)) ** 2! / legendre%coef(i)
        end do
        shape => computation_context%scattering_context%scatterer%shape
        sca = convertCtoQ(sca * 2q0 * PI / calculation_point%k**2, shape%spheroidal_type, shape%rv, shape%ab, shape%alpha%value)

    end function get_scattering_factor_far_pq

    real(knd) function kroneker(n, m)

        integer :: n, m

        if (n == m) then
            kroneker = 1q0
        else
            kroneker = 0q0
        end if
    end function kroneker

    subroutine calculate_spherical_kappa(legendre, Kappa, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Kappa(matrix_size, matrix_size)
        integer m, n, l, matrix_size

        m = legendre%m
        Kappa = 0
        !write(*,*) 'spherical kappa'
        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                Kappa(l - m + 1, n - m + 1) = -(kroneker(n + 1, l) * (l - 1) * (l - m) / (2 * l - 1) - &
                        kroneker(n - 1, l) * (l + 2) * (l + m + 1) / (2 * l + 3)) / &
                        sqrt(legendre%coef(l - m + 1) / legendre%coef(n - m + 1)) ! qsqrt((2q0 * l + 1) / (2q0 * n + 1))
                !write(*,*) 'n = ', l - m + 1, 'l = ', n - m + 1, 'val = ', Kappa(l - m + 1, n - m + 1)
            enddo
        enddo
        !write(*,*)

    end subroutine calculate_spherical_kappa

    subroutine calculate_spherical_omega(legendre, Omega, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Omega(matrix_size, matrix_size)
        integer m, n, l, matrix_size
        m = legendre%m
        Omega = 0

        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                Omega(l - m + 1, n - m + 1) = (kroneker(n, l) * 2q0 * (l * l + l + m * m - 1) / (2 * l - 1) / (2 * l + 3) - &
                        kroneker(n + 2, l) * (l - m) * (l - m - 1) / (2 * l - 1) / (2 * l - 3) - &
                        kroneker(n - 2, l) * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 5)) / &
                        sqrt(legendre%coef(l - m + 1) / legendre%coef(n - m + 1))!qsqrt((2q0 * l + 1) / (2q0 * n + 1))
            enddo
        end do
    end subroutine calculate_spherical_omega

    subroutine calculate_spherical_tau(legendre, Tau, matrix_size)
        type(LegendreCalculation) :: legendre
        complex(knd), intent(out) :: Tau(matrix_size, matrix_size)
        integer m, n, l, matrix_size

        m = legendre%m
        Tau = 0
        do n = m, m + matrix_size - 1
            do l = m, m + matrix_size - 1
                !write(*,*) 'n = ', n, 'l = ', l, 'val = ', kroneker(n, l) * l * (l + 1)
                Tau(l - m + 1, n - m + 1) = kroneker(n, l) * l * (l + 1)
            enddo
        end do
    end subroutine calculate_spherical_tau

    real(knd) function get_scattering_factor_far_uv(computation_context, mode_item, solution) result(sca)
        type(ComputationContext), intent(inout) :: computation_context
        type(ModeItem), intent(in) :: mode_item
        complex(knd), intent(in) :: solution(:)

        type(LegendreCalculation) :: legendre
        integer :: m, lnum, i, j
        type(WavelengthPoint), pointer :: calculation_point
        type(AngleType), pointer :: alpha
        type(SpheroidalShape), pointer :: shape

        complex(knd) :: ideg1, Omega(mode_item%lnum, mode_item%lnum), &
            Kappa(mode_item%lnum, mode_item%lnum), Tau(mode_item%lnum, mode_item%lnum)
        real(knd) :: k1, coef

        m = mode_item%m
        lnum = mode_item%lnum
        calculation_point => computation_context%scattering_context%calculation_point
        alpha => computation_context%scattering_context%directions%alpha
        call legendre%set(m, lnum + m - 1, alpha%angle_cos)
        call legendre%calculate()

        k1 = calculation_point%k

        sca = 0
        ideg1 = cmplx(0q0, 1q0, knd)

        call calculate_spherical_omega(legendre, Omega, lnum)
        call calculate_spherical_kappa(legendre, Kappa, lnum)
        call calculate_spherical_tau(legendre, Tau, lnum)
        do i = 1, lnum
        do j = 1, lnum
            !coef = qsqrt(legendre%coef(i) * legendre%coef(j))
            sca = sca + real(ideg1**(j - i) * (&
                    k1**2 * solution(i) * conjg(solution(j)) * Omega(i, j) + &
                            IDEG(1) * k1 * (&
                                    solution(i + lnum) * conjg(solution(j)) * &
                                            Kappa(j, i) - &
                                            conjg(solution(j + lnum)) * solution(i) *&
                                                    Kappa(i, j)) + &
                            solution(i + lnum) * conjg(solution(j + lnum)) * &
                                    Tau(i, j)), knd)! / coef
        end do
        end do
        shape => computation_context%scattering_context%scatterer%shape
        sca = convertCtoQ(sca * PI / k1 ** 2, shape%spheroidal_type, shape%rv, shape%ab, shape%alpha%value)
        if (legendre%m == 0) then
        sca = sca / 2q0
        end if

    end function get_scattering_factor_far_uv

end module spheroidal_scattering