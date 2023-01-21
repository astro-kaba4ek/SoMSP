module spheroidal_indicatrix
    use regime
    use constants
    use spheroidal
    use regime
    use spheroidal_scatterer
    use wavelength_point
    use angle
    use contexts
    implicit none

contains
    function calculate_amplitude_matrix_m(basis_type, scattering_context, m, lnum, ntheta, thetas, nphi, phis, &
        direction_calculation, solution_te, solution_tm) result(amplitude_matrix)
        integer, intent(in) :: basis_type
        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: m, lnum, ntheta, nphi
        real(knd) :: args(ntheta)
        type(AngleType), intent(in) :: thetas(ntheta), phis(nphi)
        type(AngleType) :: mphis(scattering_context%directions%nphi)
        complex(knd), intent(in) :: solution_te(:), solution_tm(:)
        complex(knd) :: amplitude_matrix(2, 2, scattering_context%directions%nphi, ntheta)

        real(knd) :: start, finish
        
        type(SpheroidalCalculation), intent(inout) :: direction_calculation
        integer :: i, j

        args = [(thetas(i)%angle_cos, i = 1, ntheta)]

        call direction_calculation%calculate(m, lnum, scattering_context%scatterer%c0(1), &
        scattering_context%scatterer%ksi(1), &
                ntheta, args, scattering_context%scatterer%spheroidal_type, .false., .true., .false.)

        mphis = [(AngleType(m * phis(i)%value), i = 1, nphi)]

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate amplitude matrix for bucket'

        100 format('#',2A8,' ',4A48)
        102 format(' ',2F8.2,' ',8E24.15)
        if (LOG_AMPLITUDE_MATRIX) then
            write(LOG_FD,*) '{INFO} calculated amplitude matrix for bucket: ', ntheta, 'thetas x', nphi, 'phis'
            write(LOG_FD,100) 'theta', &
                'phi', &
                'A11','A12','A21','A22'
        endif
        call cpu_time(start)
        do i = 1, ntheta
            do j = 1, nphi
                amplitude_matrix(:,:,j,i) = get_amplitude_matrix(&
                        basis_type, lnum, solution_te, solution_tm, scattering_context%calculation_point, &
                        thetas(i), mphis(j), direction_calculation, i)

                if (LOG_AMPLITUDE_MATRIX) then
                    write(LOG_FD,102) thetas(i)%value * 180q0 / PI, &
                        phis(j)%value * 180q0 / PI, &
                        amplitude_matrix(:,:,j,i)
                endif
            end do
        end do
        call cpu_time(finish)
        call log_time('amplitude matrix bucket', finish - start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate amplitude matrix for bucket'

    end function calculate_amplitude_matrix_m

    real(knd) function update_amplitude_and_get_accuracy(ampl, update, ntheta, nphi) result(accuracy)
        integer, intent(in) :: ntheta, nphi
        complex(knd), dimension(2,2,nphi,ntheta), intent(inout) :: ampl
        complex(knd), dimension(2,2,nphi,ntheta), intent(in) :: update
        integer :: i, j, k, n
        logical :: wrote

        ampl = ampl + update
        accuracy = 0
        wrote = .false.
        do n = 1, ntheta
            do k = 1, nphi
                do j = 1, 2
                    do i = 1, 2
                        if (abs(ampl(i, j, k, n)) > 1.0e-64_knd) then
                            accuracy = max(accuracy, abs(update(i,j,k,n) / ampl(i,j,k,n)))
                        endif
                    enddo
                enddo
            enddo
        enddo
    end function update_amplitude_and_get_accuracy

    function get_scatterring_matrix_for_bucket(scattering_context, minm, maxm, &
            matrix_size, solution_te, solution_tm, &
            ntheta, thetas, nphi, phis) result(scat)
        type(ScatteringContext), intent(in) :: scattering_context
        integer, intent(in) :: minm, maxm, matrix_size, ntheta, nphi
        complex(knd), intent(in) :: solution_te(2 * matrix_size, minm:maxm), solution_tm(2 * matrix_size, minm:maxm)
        type(AngleType), intent(in) :: thetas(ntheta), phis(nphi, minm:maxm)

        type(SpheroidalCalculation) :: direction_calculation
        real(knd) :: arg(ntheta)

        real(knd) :: scat(4,4,nphi,ntheta)
        complex(knd) :: ampl(2,2,nphi,ntheta), part(2,2)
        integer :: m, i, j

        scat = 0
        ampl = 0
        arg = (/(thetas(i)%angle_cos, i = 1, ntheta)/)
        do m = minm, maxm
            ! write(*,*) 'allocating for direction calculation: ms = ', matrix_size, ' nargs = ', size(arg)
            call direction_calculation%calculate(m, matrix_size, scattering_context%scatterer%c0(1), &
            scattering_context%scatterer%ksi(1), &
                    ntheta, arg, scattering_context%scatterer%spheroidal_type, .false., .true., .false.)
            do i = 1, ntheta
                do j = 1, nphi
                    part = 0
                    part = get_amplitude_matrix_uv(&
                            matrix_size, solution_te(:,m), solution_tm(:,m), scattering_context%calculation_point, &
                            thetas(i), phis(j, m), direction_calculation, i)
                    ampl(:,:,j,i) = ampl(:,:,j,i) + part
                    ! write(*,*) 'm = ', m, ' part = ', part, ' dpart = ', part / ampl(:,:,j,i)
                end do
            end do
        end do
        do i = 1, ntheta
            do j = 1, nphi
                scat(:,:,j,i) = get_scattering_matrix(ampl(:,:,j,i))
            end do
        end do
        !write(*,*) 'scat = ', scat
    end function get_scatterring_matrix_for_bucket

    complex(knd) function get_axisymmetric_amplitude(matrix_size, solution, direction_calculation, argnum) result(T11)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(matrix_size)

        integer :: i
        !write(*, *) 'get_axisymmetric_amplitude_1_1', ' m = ', this%m
        T11 = 0

        do i = 1, matrix_size
            T11 = T11 + NEGIDEG(mod(i, 4)) * solution(i) * direction_calculation%s1(i, argnum)
        end do

    end function get_axisymmetric_amplitude

    complex(knd) function get_nonaxisymmetric_amplitude_1_1(matrix_size, solution, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(T11)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i

        !write(*, *) 'get_nonaxisymmetric_amplitude_1_1', ' m = ', this%m

        T11 = 0

        do i = 1, matrix_size
            T11 = T11 + NEGIDEG(mod(i + direction_calculation%m + 2, 4)) * &
                    (calculation_point%k * solution(i) * direction_calculation%s1(i, argnum) + &
                            IDEG(1) * solution(i + matrix_size) * direction_calculation%s1d(i, argnum))
        end do

        T11 = T11 * theta%angle_sin * mphi%angle_cos

    end function get_nonaxisymmetric_amplitude_1_1

    complex(knd) function get_nonaxisymmetric_amplitude_1_2(matrix_size, solution, &
            theta, mphi, direction_calculation, argnum) result(T12)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        integer :: i
        !write(*,*) 'get_nonaxisymmetric_amplitude_1_2', ' m = ', this%m

        T12 = 0
        do i = 1, matrix_size
            T12 = T12 - NEGIDEG(mod(i + direction_calculation%m + 3, 4)) * solution(i + matrix_size) * &
                    direction_calculation%m * direction_calculation%s1(i, argnum)
        end do

        T12 = T12 / theta%angle_sin * mphi%angle_sin

    end function get_nonaxisymmetric_amplitude_1_2

    complex(knd) function get_nonaxisymmetric_amplitude_2_1(matrix_size, solution, &
            theta, mphi, direction_calculation, argnum) result(T21)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)

        !write(*, *) 'get_nonaxisymmetric_amplitude_2_1', ' m = ', this%m
        T21 = get_nonaxisymmetric_amplitude_1_2(matrix_size, solution, &
                theta, mphi, direction_calculation, argnum)

    end function get_nonaxisymmetric_amplitude_2_1

    complex(knd) function get_nonaxisymmetric_amplitude_2_2(matrix_size, solution, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(T22)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution(2 * matrix_size)
        !write(*, *) 'get_nonaxisymmetric_amplitude_2_2', ' m = ', this%m

        T22 = -get_nonaxisymmetric_amplitude_1_1(matrix_size, solution, calculation_point, &
                theta, mphi, direction_calculation, argnum)

    end function get_nonaxisymmetric_amplitude_2_2

    function get_amplitude_matrix(basis_type, lnum, solution_te, solution_tm, calculation_point, &
        theta, mphi, direction_calculation, argnum) result(ampl)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: lnum, argnum, basis_type
        complex(knd), intent(in) :: solution_te(:), solution_tm(:)

        complex(knd) :: ampl(2,2)

        if (basis_type == UV) then
            ampl = get_amplitude_matrix_uv(lnum, solution_te, solution_tm, calculation_point, &
            theta, mphi, direction_calculation, argnum)
        elseif (basis_type == PQ) then
            ampl = get_amplitude_matrix_pq(lnum, solution_te, solution_tm, direction_calculation, argnum)
        else
            call assert(.false., 'unknown basis type for amplitude matrix')
        endif

    end function get_amplitude_matrix

    function get_amplitude_matrix_uv(matrix_size, solution_te, solution_tm, calculation_point, &
            theta, mphi, direction_calculation, argnum) result(ampl)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        type(WavelengthPoint), intent(in) :: calculation_point
        type(AngleType), intent(in) :: theta, mphi
        integer, intent(in) :: matrix_size, argnum
        complex(knd), intent(in) :: solution_te(2 * matrix_size), solution_tm(2 * matrix_size)

        complex(knd) :: ampl(2,2)

        ampl(1,1) = get_nonaxisymmetric_amplitude_1_1(matrix_size, solution_tm, calculation_point, &
                theta, mphi, direction_calculation, argnum)
        ampl(1,2)= get_nonaxisymmetric_amplitude_1_2(matrix_size, solution_te, &
                theta, mphi, direction_calculation, argnum)
        ampl(2,1) = get_nonaxisymmetric_amplitude_2_1(matrix_size, solution_tm, &
                theta, mphi, direction_calculation, argnum)
        ampl(2,2) = get_nonaxisymmetric_amplitude_2_2(matrix_size, solution_te, calculation_point, &
                theta, mphi, direction_calculation, argnum)
        if (direction_calculation%m == 0) then
            ampl = ampl * 0.5_knd
        end if
        !write(*,*) 'theta = ', theta%value, 'phi = ', mphi%value, 'ampl = ', ampl
        ! write(*,*) 'uv ampl = ', ampl

    end function get_amplitude_matrix_uv

    function get_amplitude_matrix_pq(lnum, solution_te, solution_tm, direction_calculation, argnum) result(ampl)
        type(SpheroidalCalculation), intent(in) :: direction_calculation
        integer, intent(in) :: lnum, argnum
        complex(knd), intent(in) :: solution_te(lnum), solution_tm(lnum)

        complex(knd) :: ampl(2,2)

        ampl(1,1) = -get_axisymmetric_amplitude(lnum, solution_tm, direction_calculation, argnum)
        ampl(1,2) = 0
        ampl(2,1) = 0
        ampl(2,2) = get_axisymmetric_amplitude(lnum, solution_te, direction_calculation, argnum)

        ! write(*,*) 'pq ampl = ', ampl

        ! write(*,*) 'theta = ', theta%value, 'phi = ', mphi%value, 'ampl = ', ampl

    end function get_amplitude_matrix_pq

    function get_scattering_matrix(ampl) result(scat)
        complex(knd), intent(in) :: ampl(2,2)

        real(knd) :: scat(4,4)

        ! write(*,*) 'ampl0 = ', ampl
        scat(1,1) = (abs(ampl(1,1))**2 + abs(ampl(1,2))**2 + abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(1,2) = (-abs(ampl(1,1))**2 - abs(ampl(1,2))**2 + abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(1,3) = real(ampl(2,2) * conjg(ampl(1,2)) + ampl(1,1) * conjg(ampl(2,1)),knd)
        scat(1,4) = imag(ampl(2,2) * conjg(ampl(1,2)) - ampl(1,1) * conjg(ampl(2,1)))
        scat(2,1) = (-abs(ampl(1,1))**2 + abs(ampl(1,2))**2 - abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(2,2) = (abs(ampl(1,1))**2 - abs(ampl(1,2))**2 - abs(ampl(2,1))**2 + abs(ampl(2,2))**2) * 0.5q0
        scat(2,3) = real(ampl(2,2) * conjg(ampl(1,2)) - ampl(1,1) * conjg(ampl(2,1)),knd)
        scat(2,4) = imag(ampl(2,2) * conjg(ampl(1,2)) + ampl(1,1) * conjg(ampl(2,1)))
        scat(3,1) = real(ampl(2,2) * conjg(ampl(2,1)) + ampl(1,1) * conjg(ampl(1,2)),knd)
        scat(3,2) = real(ampl(2,2) * conjg(ampl(2,1)) - ampl(1,1) * conjg(ampl(1,2)),knd)
        scat(3,3) = real(ampl(1,1) * conjg(ampl(2,2)) + ampl(1,2) * conjg(ampl(2,1)),knd)
        scat(3,4) = imag(ampl(2,2) * conjg(ampl(1,1)) + ampl(2,1) * conjg(ampl(1,2)))
        scat(4,1) = imag(ampl(2,1) * conjg(ampl(2,2)) + ampl(1,1) * conjg(ampl(1,2)))
        scat(4,2) = imag(ampl(2,1) * conjg(ampl(2,2)) - ampl(1,1) * conjg(ampl(1,2)))
        scat(4,3) = imag(ampl(1,1) * conjg(ampl(2,2)) - ampl(1,2) * conjg(ampl(2,1)))
        scat(4,4) = real(ampl(1,1) * conjg(ampl(2,2)) - ampl(1,2) * conjg(ampl(2,1)),knd)
        ! write(*,*) 'ampl = ', ampl
        ! write(*,*) 'scat = ', scat

    end function get_scattering_matrix

end module spheroidal_indicatrix