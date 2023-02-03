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

end module spheroidal_scattering