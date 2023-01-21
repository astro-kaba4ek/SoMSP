! Created by drakosha on 13.07.2021.

module spheroidal_tmatrix
    use regime
    use wavelength_point
    use spheroidal_scatterer
    use integrals
    use spheroidal
    use matrix
    use contexts
    use logging
!    use spherical
    implicit none

contains
    ! formulas (23)-(26) AppliedOptics 2012
    ! A_ik = W(left_R Delta - mu_j/mu_{j+1}Delta right_R - (mu_j/mu_{j+1}-1)ksi/(ksi^20f)Delta)P
    function get_part(f, ksi, mu, left_R, right_R, W, P, &
            Delta, matrix_size) result(res)
        complex(knd) :: mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: left_R(matrix_size), right_R(matrix_size), W(matrix_size), P(matrix_size, matrix_size), &
                Delta(matrix_size, matrix_size), tmp(matrix_size, matrix_size), res(matrix_size, matrix_size), val

        tmp = Delta
        call multiply_by_diag_left(tmp, matrix_size, left_R)
        res = tmp

        tmp = Delta
        call multiply_by_diag_right(tmp, matrix_size, right_R)
        res = res - mu(0) / mu(1) * tmp - (mu(0) / mu(1) - 1.0_knd)*ksi / (ksi**2 - f) * Delta
        call multiply_by_diag_left(res, matrix_size, W)

        res = matmul(res, P)
    end function get_part

    subroutine calculate_tmatrix_sph_pq(scatterering_context, computation_context, &
            is_te_mode, tmatrix)

        type(ScatteringContext), intent(in) :: scatterering_context
        type(SpheroidalContext), intent(in) :: computation_context
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: tmatrix(computation_context%lnum, computation_context%lnum)

        integer :: lnum, nol
        integer :: n, i, l, j
        ! by layer
        complex(knd) :: mu(0:computation_context%nol)
        ! by expansion
        complex(knd), dimension(computation_context%lnum) :: R11, R31, R12, W1, R32
        complex(knd), dimension(computation_context%lnum, computation_context%lnum) :: adder, A11, A31, A31inv, P
        complex(knd) :: big_matr(2 * computation_context%lnum, computation_context%lnum)
        complex(knd) :: tmp(2 * computation_context%lnum, 2 * computation_context%lnum)

        real(knd), allocatable :: mult_coef(:)

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate tmatrix sph pq'

        nol = computation_context%nol
        lnum = computation_context%lnum

        mu = scatterering_context%calculation_point%mu
        if (.not. is_te_mode) then
            mu = scatterering_context%calculation_point%eps
        endif

        ! setting up core
        R11 = computation_context%layers(0, nol)%r1d(1:lnum) / computation_context%layers(0, nol)%r1(1:lnum)
        R31 = computation_context%layers(0, nol)%r3d(1:lnum) / computation_context%layers(0, nol)%r3(1:lnum)
        R12 = computation_context%layers(1, nol)%r1d(1:lnum) / computation_context%layers(1, nol)%r1(1:lnum)
        W1 = -1.0_knd / (R31 - R11)
        call get_identity_matrix(P, lnum)

        big_matr = 0
        big_matr(lnum + 1:,:) = get_part( &
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(nol), &
                mu(nol - 1:nol), &
                R11, R12, W1, &
                computation_context%Pi1(:,:,nol), &
                computation_context%Delta(:,:,nol), &
                lnum)
        big_matr(1 : lnum,:) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(nol), &
                mu(nol - 1:nol), &
                R31, R12, W1, &
                computation_context%Pi1(:,:,nol), &
                computation_context%Delta(:,:,nol), &
                lnum)

        do j = nol - 1, 1, -1
            R11 = computation_context%layers(0, j)%r1d(1:lnum) / computation_context%layers(0, j)%r1(1:lnum)
            R31 = computation_context%layers(0, j)%r3d(1:lnum) / computation_context%layers(0, j)%r3(1:lnum)
            R12 = computation_context%layers(1, j)%r1d(1:lnum) / computation_context%layers(1, j)%r1(1:lnum)
            R32 = computation_context%layers(1, j)%r3d(1:lnum) / computation_context%layers(1, j)%r3(1:lnum)
            W1 = -1.0_knd / (R31 - R11)
            tmp = 0

            tmp(lnum + 1: 2*lnum, 1:lnum) = get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R11, R12, W1, &
                computation_context%Pi1(:,:,j), &
                computation_context%Delta(:,:,j), &
                lnum)
            tmp(1:lnum, 1:lnum) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R31, R12, W1, &
                computation_context%Pi1(:,:,j), &
                computation_context%Delta(:,:,j), &
                lnum)

            tmp(lnum + 1: 2 * lnum, lnum + 1:2 * lnum) = get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R11, R32, W1, &
                computation_context%Pi3(:,:,j), &
                computation_context%Delta(:,:,j), &
                lnum)
            tmp(1:lnum, lnum + 1:) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R31, R32, W1, &
                computation_context%Pi3(:,:,j), &
                computation_context%Delta(:,:,j), &
                lnum)
            big_matr = matmul(tmp, big_matr)
        enddo

        A11 = 0
        A31 = 0
        A31inv = 0
        A31 = big_matr(1:lnum, :)
        A11 = big_matr(lnum + 1:2 * lnum, :)

        call multiply_by_diag_left(A11, lnum, 1.0_knd / computation_context%layers(0,1)%r3(1:lnum))
        call quick_inverse_matrix(A31, lnum, A31inv)
        call multiply_by_diag_right(A31inv, lnum, computation_context%layers(0,1)%r1(1:lnum))

        tmatrix = -matmul(A11, A31inv)

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate tmatrix sph pq'

    end subroutine calculate_tmatrix_sph_pq

    !  Nonaxisymmetric part

    function get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Epsilon, matrix_size) result(res)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                res(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        res = -(eps(1) / eps(0) - mu(0) / mu(1)) * f * ksi / (ksi ** 2 - f) * matmul(Q01, Epsilon) - &
                (eps(1) / eps(0) - 1q0) * ksi * (Q01 - 2q0 * ksi**2 * Q01Q11)

        adder = (mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - mu(0) / mu(1) * Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        res = res + adder

        adder = Delta + (eps(1) / eps(0) - 1q0) * ksi**2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        res = res + adder

    end function get_part_11

    function get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1), identity(matrix_size, matrix_size)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = -(eps(1) / eps(0) - mu(0) / mu(1)) * f / (ksi ** 2 - f) * &
                (matmul((ksi**2 * Q01 - Delta), Kappa) + matmul(Delta, Gamma11)) + &
                (eps(1) / eps(0) - 1q0) * f * ksi**2 * 2q0 * matmul(Q01Q11, Gamma11)
        adder = (mu(0) / mu(1) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = (eps(1) / eps(0) - 1q0) * f * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_12

    function get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi**2 / (ksi ** 2 - f) * matmul(Q01, Kappa) - &
                (eps(1) / eps(0) - 1q0) * ksi**2 * 2q0 * matmul(Q01Q11, Gamma11)
        adder = -(mu(0) / mu(1) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = -(eps(1) / eps(0) - 1q0) * ksi * matmul(Q01, Gamma11)
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_21

    function get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            Delta, Q01, Q11, Q01Q11, Epsilon, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(matrix_size), right_multiplier(matrix_size), &
                adder(matrix_size, matrix_size), Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size), identity(matrix_size, matrix_size), &
                result(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size)

        call get_identity_matrix(identity, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi / (ksi ** 2 - f) * (f * matmul(Q01, Epsilon) + Delta) + &
                (eps(1) / eps(0) - 1q0) * ksi * (Q01 - 2q0 * ksi**2 * Q01Q11)

        adder = -(mu(0) / mu(1) - 1q0) * ksi**2 * Q01 - Delta
        call multiply_by_diag_right(adder, matrix_size, right_multiplier)
        result = result + adder

        adder = eps(1) / eps(0) * Delta - (eps(1) / eps(0) - 1q0) * ksi**2 * Q01
        call multiply_by_diag_left(adder, matrix_size, left_multiplier)
        result = result + adder

    end function get_part_22

    subroutine set_full_matrix(f, ksi, eps, mu, first_multiplier, left_multiplier, right_multiplier, last_multiplier, &
            Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, matrix_size, result, ddd)

        complex(knd), intent(in) :: eps(0:1), mu(0:1)
        real(knd), intent(in) :: ksi
        integer, intent(in) :: matrix_size, f
        real(knd), intent(in), optional :: ddd
        complex(knd), intent(in) :: first_multiplier(matrix_size), left_multiplier(matrix_size), right_multiplier(matrix_size), &
                Delta(matrix_size, matrix_size), Q01(matrix_size, matrix_size), &
                Q11(matrix_size, matrix_size), Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), &
                Epsilon(matrix_size, matrix_size), &
                Q01Q11(matrix_size, matrix_size), last_multiplier(matrix_size, matrix_size)
        complex(knd), intent(out) :: result(2 * matrix_size, 2 * matrix_size)
        
        integer :: i, j


        result = 0
        result(1:matrix_size, 1:matrix_size) = &
                get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Epsilon, matrix_size)
        result(1:matrix_size, (matrix_size + 1):(2 * matrix_size)) = &
                get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), 1:matrix_size) = &
                get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Kappa, Gamma11, matrix_size)
        result((matrix_size + 1):(2 * matrix_size), (matrix_size + 1):(2 * matrix_size)) = &
                get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, Delta, Q01, Q11, Q01Q11, &
                        Epsilon, matrix_size)
        do i = 0, 1
            do j = 0, 1
                call multiply_by_diag_left(&
                        result((i * matrix_size + 1):((i + 1) * matrix_size), (j * matrix_size + 1):((j + 1) * matrix_size)), &
                        matrix_size, first_multiplier)
                result((i * matrix_size + 1):((i + 1) * matrix_size), (j * matrix_size + 1):((j + 1) * matrix_size)) =&
                        matmul(result((i * matrix_size + 1):((i + 1) * matrix_size), &
                                (j * matrix_size + 1):((j + 1) * matrix_size)), last_multiplier)
            end do
        end do
        if (present(ddd)) then
            result(:, matrix_size + 1 : 2 * matrix_size) = ddd * result(:, matrix_size + 1 : 2 * matrix_size)
        end if

    end subroutine set_full_matrix

    subroutine calculate_tmatrix_spheroidal_uv(scattering_context, computation_context, is_te_mode, tmatrix)

        type(ScatteringContext), intent(in) :: scattering_context
        type(SpheroidalContext), intent(in) :: computation_context
        logical, intent(in) :: is_te_mode
        complex(knd), intent(out) :: tmatrix(2 * computation_context%lnum, 2 * computation_context%lnum)

        integer :: n, i, l, nol, j, lnum
        real(knd) :: k1, ddd, start, finish, global_start
        complex(knd) :: c1
        complex(knd), dimension(0:scattering_context%scatterer%number_of_layers)  :: mu, eps
        complex(knd), dimension(computation_context%lnum) :: R11, R31, R12, R32, W1
        complex(knd), dimension(2 * computation_context%lnum) :: initial_corrector, solution_corrector
        complex(knd), dimension(2 * computation_context%lnum, 2 * computation_context%lnum) :: A11, A31, A31inv
        complex(knd), dimension(4 * computation_context%lnum, 2 * computation_context%lnum) :: big_matr
        complex(knd), dimension(4 * computation_context%lnum, 4 * computation_context%lnum) :: tmp

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate tmatrix sph uv'

        call cpu_time(global_start)
        mu = scattering_context%calculation_point%mu
        eps = scattering_context%calculation_point%eps
        if (.not. is_te_mode) then
            mu = scattering_context%calculation_point%eps
            eps = scattering_context%calculation_point%mu
        endif

        k1 = scattering_context%calculation_point%k
        c1 = scattering_context%scatterer%c0(1)
        nol = scattering_context%scatterer%number_of_layers
        lnum = computation_context%lnum

        ! setting up core
        R11 = computation_context%layers(0, nol)%r1d(1:lnum) / computation_context%layers(0, nol)%r1(1:lnum)
        R31 = computation_context%layers(0, nol)%r3d(1:lnum) / computation_context%layers(0, nol)%r3(1:lnum)
        R12 = computation_context%layers(1, nol)%r1d(1:lnum) / computation_context%layers(1, nol)%r1(1:lnum)
        W1 = -1.0_knd / (R31 - R11)

        big_matr = 0
        call cpu_time(start)
        ! A31
        call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(nol), &
                eps(nol - 1: nol), &
                mu(nol - 1: nol), &
                W1, &
                R31, &
                R12, &
                computation_context%Pi1(:,:,nol), &
                computation_context%Delta(:,:,nol), &
                computation_context%Q01(:,:,nol), &
                computation_context%Q11(:,:,nol), &
                computation_context%Q01Q11(:,:,nol), &
                computation_context%Kappa(:,:,nol), &
                computation_context%Gamma11(:,:,nol), &
                computation_context%Epsilon(:,:,nol), &
                lnum, &
                big_matr(1 : 2 * lnum,:))
        ! A11
        call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(nol), &
                eps(nol - 1: nol), &
                mu(nol - 1: nol), &
                W1, &
                R11, &
                R12, &
                computation_context%Pi1(:,:,nol), &
                computation_context%Delta(:,:,nol), &
                computation_context%Q01(:,:,nol), &
                computation_context%Q11(:,:,nol), &
                computation_context%Q01Q11(:,:,nol), &
                computation_context%Kappa(:,:,nol), &
                computation_context%Gamma11(:,:,nol), &
                computation_context%Epsilon(:,:,nol), &
                lnum, &
                big_matr(2 * lnum + 1 : 4 * lnum, :))
        big_matr(1 : 2 * lnum,:) = -big_matr(1 : 2 * lnum,:)
        call cpu_time(finish)
        call log_time('tmatrix core', finish - start)    

        call cpu_time(start)
        do j = nol - 1, 1, -1
            R11 = computation_context%layers(0, j)%r1d(1:lnum) / computation_context%layers(0, j)%r1(1:lnum)
            R31 = computation_context%layers(0, j)%r3d(1:lnum) / computation_context%layers(0, j)%r3(1:lnum)
            R12 = computation_context%layers(1, j)%r1d(1:lnum) / computation_context%layers(1, j)%r1(1:lnum)
            R32 = computation_context%layers(1, j)%r3d(1:lnum) / computation_context%layers(1, j)%r3(1:lnum)
            W1 = -1.0_knd / (R31 - R11)
            ddd = scattering_context%scatterer%d(j) / scattering_context%scatterer%d(j + 1)

            tmp = 0

            ! A31
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R31, &
                R12, &
                computation_context%Pi1(:,:,j), &
                computation_context%Delta(:,:,j), &
                computation_context%Q01(:,:,j), &
                computation_context%Q11(:,:,j), &
                computation_context%Q01Q11(:,:,j), &
                computation_context%Kappa(:,:,j), &
                computation_context%Gamma11(:,:,j), &
                computation_context%Epsilon(:,:,j), &
                lnum, &
                tmp(1:2 * lnum, 1:2 * lnum), &
                ddd)

            ! A11
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R11, &
                R12, &
                computation_context%Pi1(:,:,j), &
                computation_context%Delta(:,:,j), &
                computation_context%Q01(:,:,j), &
                computation_context%Q11(:,:,j), &
                computation_context%Q01Q11(:,:,j), &
                computation_context%Kappa(:,:,j), &
                computation_context%Gamma11(:,:,j), &
                computation_context%Epsilon(:,:,j), &
                lnum, &
                tmp(2 * lnum + 1: 4 * lnum, 1:2 * lnum), &
                ddd)

            ! A33
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R31, &
                R32, &
                computation_context%Pi3(:,:,j), &
                computation_context%Delta(:,:,j), &
                computation_context%Q01(:,:,j), &
                computation_context%Q11(:,:,j), &
                computation_context%Q01Q11(:,:,j), &
                computation_context%Kappa(:,:,j), &
                computation_context%Gamma11(:,:,j), &
                computation_context%Epsilon(:,:,j), &
                lnum, &
                tmp(1:2 * lnum, 2 * lnum + 1: 4 * lnum), &
                ddd)
            ! A13
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R11, &
                R32, &
                computation_context%Pi3(:,:,j), &
                computation_context%Delta(:,:,j), &
                computation_context%Q01(:,:,j), &
                computation_context%Q11(:,:,j), &
                computation_context%Q01Q11(:,:,j), &
                computation_context%Kappa(:,:,j), &
                computation_context%Gamma11(:,:,j), &
                computation_context%Epsilon(:,:,j), &
                lnum, &
                tmp(2 * lnum + 1:4 * lnum, 2 * lnum + 1: 4 * lnum), &
                ddd)
            tmp(1 : 2 * lnum,:) = -tmp(1 : 2 * lnum,:)
            big_matr = matmul(tmp, big_matr)

        enddo
        call cpu_time(finish)
        call log_time('tmatrix layers', finish - start)

        call cpu_time(start)
        A11 = 0
        A31 = 0
        A31inv = 0
        A31 = big_matr(1:2 * lnum, :)
        A11 = big_matr(2 * lnum + 1:4 * lnum, :)

        initial_corrector(1:lnum) = 1q0 / (computation_context%layers(0, 1)%r1(1:lnum) * k1)
        initial_corrector((lnum + 1):(2 * lnum)) = 1q0 / (computation_context%layers(0, 1)%r1(1:lnum) * c1)
        solution_corrector(1:lnum) = 1q0 / (computation_context%layers(0, 1)%r3(1:lnum) * k1)
        solution_corrector((lnum + 1):(2 * lnum)) = 1q0 / (computation_context%layers(0, 1)%r3(1:lnum) * c1)

        call multiply_by_diag_left(A31, 2 * lnum, initial_corrector)
        call multiply_by_diag_left(A11, 2 * lnum, solution_corrector)
        call cpu_time(finish)
        call log_time('tmatrix basis correction', finish - start)
        

        call cpu_time(start)
        call quick_inverse_matrix(A31, 2 * lnum, A31inv)
        call cpu_time(finish)
        call log_time('tmatrix inversion', finish - start)

        call cpu_time(start)
        tmatrix = matmul(A11, A31inv)
        call cpu_time(finish)
        call log_time('tmatrix multiplication', finish - start)
        call log_time('tmatrix total', finish - global_start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate tmatrix sph uv'

    end subroutine calculate_tmatrix_spheroidal_uv


end module spheroidal_tmatrix