! Created by odawing on 13.08.22.

module integrals
    use regime
    use spheroidal
    use constants
    implicit none
    integer, parameter :: bucket = 20
    real(knd), parameter :: accuracy = 1d-32
    real(knd), parameter :: border = 1.0e-128_knd

    abstract interface
        real(knd) pure function func_coef(m, n)
            import :: knd
            integer, intent(in) :: m, n
        end function func_coef
    end interface
contains
    real(knd) pure function delta_coef(m, n)
        integer, intent(in) :: m, n
        delta_coef = 1.0_knd
    end function delta_coef

    logical function update_with_accuracy(value, delta, accuracy)
        complex(knd), intent(inout) :: value
        complex(knd), intent(in) :: delta
        real(knd), intent(in) :: accuracy

        value = value + delta
        update_with_accuracy = (abs(delta / value) < accuracy)
    end function update_with_accuracy

    !  fills the array common_multiplier with real numbers
    !  (r+2m)/r!*2/(2r+2m+1) for r = 0..accuracy
    !  common_multiplier is allocated if necessary
    subroutine fill_common_multiplier(m, accuracy, common_multiplier)
        integer :: m, accuracy, r
        real(knd), allocatable, dimension(:) :: common_multiplier
        real(knd) :: factorial

        if (allocated(common_multiplier) .and. size(common_multiplier) < accuracy + 1) then
            deallocate(common_multiplier)
        end if
        if (.not. allocated(common_multiplier)) then
            allocate(common_multiplier(0:accuracy))
        end if

        factorial = 1q0
        do r = 2, 2 * m
            factorial = factorial * r
        end do
        do r = 0, accuracy
            common_multiplier(r) = factorial * 2q0 / (2 * r + 2 * m + 1)
            factorial = factorial * (r + 2 * m + 1) / (r + 1)
        enddo
    end subroutine fill_common_multiplier

    !  Delta_{n,l}^{(m)}(c_i, c_j) = \int_{-1}^{1} S_{m,n}(c_i, \eta) S_{m,n}(c_j, \eta) d\eta =
    !  \sum'_{r=0}^{accuracy} d_r^{m,n}(c_i) d_r^{m,l}(c_j) common_multiplier(m, r)
    !  in sum' only r = n-l (mod 2) are used,
    !  angular spheroidal functions and legendre coefficients are normalized
    !  and taken from SpheroidalCalculation variables
    !  common_multiplier is calculated in subroutine fill_common_multiplier
    real(knd) pure function tau_c_middle(m, r)
        integer, intent(in) :: m, r
        tau_c_middle = real((r + m) * (r + m + 1), knd)
    end function tau_c_middle

    pure real(knd) function gamma_c_lower(m, r)
        integer, intent(in) :: m, r
        gamma_c_lower = real(r, knd) / (2 * r + 2 * m - 1)
    end function gamma_c_lower
    pure real(knd) function gamma_c_upper(m, r)
        integer, intent(in) :: m, r
        gamma_c_upper = real(2 * m + 1 + r, knd) / (2 * r + 2 * m + 3)
    end function gamma_c_upper

    pure real(knd) function kappa_c_lower(m, r)
        integer, intent(in) :: m, r
        kappa_c_lower = -real(r * (r + m - 1), knd) / (2 * r + 2 * m - 1)
    end function kappa_c_lower
    pure real(knd) function kappa_c_upper(m, r)
        integer, intent(in) :: m, r
        kappa_c_upper = real((r + m + 2) * (2 * m + 1 + r), knd) / (2 * r + 2 * m + 3)
    end function kappa_c_upper

    ! Sigma
    pure real(knd) function sigma_c_lower(m, r)
        integer, intent(in) :: m, r
        sigma_c_lower = -real((r - 1) * r * (r + m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function sigma_c_lower
    pure real(knd) function sigma_c_middle(m, r)
        integer, intent(in) :: m, r
        sigma_c_middle = real(3 * (r + m) * (r + m + 1) - m * m - 2, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function sigma_c_middle
    pure real(knd) function sigma_c_upper(m, r)
        integer, intent(in) :: m, r
        sigma_c_upper = real((r + m + 2) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function sigma_c_upper

    pure real(knd) function epsilon_c_lower(m, r)
        integer, intent(in) :: m, r
        epsilon_c_lower = -real((r - 1) * r * (r + m - 2), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function epsilon_c_lower
    pure real(knd) function epsilon_c_middle(m, r)
        integer, intent(in) :: m, r
        epsilon_c_middle = real((r + m) * (r + m + 1) - 3 * m * m, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function epsilon_c_middle
    pure real(knd) function epsilon_c_upper(m, r)
        integer, intent(in) :: m, r
        epsilon_c_upper = real((r + m + 3) * (r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function epsilon_c_upper

    ! Omega
    pure real(knd) function omega_c_lower(m, r)
        integer, intent(in) :: m, r
        omega_c_lower = -real((r - 1) * r, knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m - 3))
    end function omega_c_lower
    pure real(knd) function omega_c_middle(m, r)
        integer, intent(in) :: m, r
        omega_c_middle = real(2 * ((r + m) * (r + m + 1) + m * m - 1), knd) / &
                ((2 * r + 2 * m - 1) * (2 * r + 2 * m + 3))
    end function omega_c_middle
    pure real(knd) function omega_c_upper(m, r)
        integer, intent(in) :: m, r
        omega_c_upper = -real((r + 2 * m + 1) * (r + 2 * m + 2), knd) / &
                ((2 * r + 2 * m + 3) * (2 * r + 2 * m + 5))
    end function omega_c_upper

    subroutine single_dep_integral(Int, m, First, Second, mult_coef, middle)
        complex(knd), intent(out) :: Int(:,:)
        integer, intent(in) :: m
        complex(knd), intent(in) :: First(0:,:), Second(0:,:)
        real(knd), intent(in) :: mult_coef(0:)
        real(knd) :: fc(0:10000)
        procedure(func_coef) :: middle

        integer n, l, r, k, ms, md
        real :: start, finish
        complex(knd) :: value

        Int = 0
        ms = size(Int(:, 1))
        md = min(size(First(:, 1)), size(Second(:, 1))) - 1
        fc = 0
        do n = 0, md
            fc(n) =  mult_coef(n) * middle(m, n)
        end do

        !write(*,*) 'common mult = ', common_multiplier(1:100)
        call cpu_time(start)
        do l = 1, ms
            do n = 2 - mod(l, 2), ms, 2
                !  sum from mod(n + 1, 2) because here n = n0 - m + 1 =>
                !  mod(n0-m, 2) = mod(n0 - m + 1 + 1, 2) = mod(n + 1, 2)
                do r = mod(n + 1, 2), md - 2 * bucket + 1, 2 * bucket
                    do k = r, r + 2 * (bucket - 1), 2
                        value = First(k, n) * Second(k, l) * fc(k)
!                        write(*,*) 'value = ', value, 'first = ', first(k, n), 'second = ', second(k, l),&
!                                'mult_coef = ', mult_coef(k), 'func_coef = ', middle(m, k)
!                                'same_assoc_value = ', First(k, n) * Second(k, l) * mult_coef(k) * middle(m, k),&
!                                'diff_assoc_value = ', First(k, n) * Second(k, l) * (mult_coef(k) * middle(m, k))
                        Int(n, l) = Int(n, l) + value
                    end do
                    if ((abs(value) < accuracy * abs(Int(n, l))).or.(abs(value) < border * abs(Int(1,1)))) then
                        !                         write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                        exit
                    end if
                enddo
            enddo
        enddo
        call cpu_time(finish)
        ! write(*,*) 'time = ', finish - start
        !        call log_matrix(FILE_DESCRIPTOR(WARNING), 'Delta', Delta, .false., matrix_size)
        !write(*,*)

    end subroutine single_dep_integral

    subroutine double_dep_integral(Int, m, First, Second, mult_coef, lower, upper)
        complex(knd), intent(out) :: Int(:,:)
        integer, intent(in) :: m
        complex(knd), intent(in) :: First(0:,:), Second(0:,:)
        real(knd), intent(in) :: mult_coef(0:)
        real(knd) :: fcoef_lower(0:20000), fcoef_upper(0:20000)
        procedure(func_coef) :: lower, upper

        integer n, l, r, k, ms, md, ix
        real :: start, finish
        complex(knd) :: value

        Int = 0
        ms = size(Int(:, 1))
        md = min(size(First(:, 1)), size(Second(:, 1))) - 1
        ! write(*,*) 'here md = ', md
        do n = 0, md
            ! write(*,*) n, size(mult_coef), size(fcoef_lower)
            fcoef_lower(n) = lower(m, n) * mult_coef(n)
            fcoef_upper(n) = upper(m, n) * mult_coef(n)
        end do

        !write(*,*) 'common mult = ', common_multiplier(1:100)
        call cpu_time(start)
        do l = 1, ms
            do n = 1 + mod(l, 2), ms, 2
                ix = mod(n + 1, 2)
                if (ix == 0) then
                    Int(n, l) = First(ix, n) * Second(ix + 1, l) * fcoef_upper(ix)
                    ix = 2
                endif

                do r = ix, md - 2 * bucket + 1, 2 * bucket
                    do k = r, r + 2 * (bucket - 1), 2
                        value = First(k, n) * (Second(k - 1, l) * fcoef_lower(k) + Second(k + 1, l) * fcoef_upper(k))
                        !                        write(*,*) 'value = ', value, 'first = ', first(k, n), 'second = ', second(k, l),&
                        !                                'mult_coef = ', mult_coef(k), 'func_coef = ', func_coef(m, k)
                        Int(n, l) = Int(n, l) + value
                    end do
                    if ((abs(value) < accuracy * abs(Int(n, l))).or.(abs(value) < border * abs(Int(1,1)))) then
                        !                         write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                        exit
                    end if
                enddo
            enddo
        enddo
        call cpu_time(finish)
        ! write(*,*) 'time = ', finish - start
        !        call log_matrix(FILE_DESCRIPTOR(WARNING), 'Delta', Delta, .false., matrix_size)
        !write(*,*)

    end subroutine double_dep_integral

    subroutine triple_dep_integral(Int, m, First, Second, mult_coef, lower, middle, upper)
        complex(knd), intent(out) :: Int(:,:)
        integer, intent(in) :: m
        complex(knd), intent(in) :: First(0:,:), Second(0:,:)
        real(knd), intent(in) :: mult_coef(0:)
        real(knd) :: fcoef_lower(0:10000), fcoef_middle(0:100000), fcoef_upper(0:10000)
        procedure(func_coef) :: lower, middle, upper

        integer n, l, r, k, ms, md, ix
        real :: start, finish
        complex(knd) :: value

        Int = 0
        ms = size(Int(:, 1))
        md = min(size(First(:, 1)), size(Second(:, 1))) - 1
        do n = 0, md
            fcoef_lower(n) = lower(m, n) * mult_coef(n)
            fcoef_middle(n) = middle(m, n) * mult_coef(n)
            fcoef_upper(n) = upper(m, n) * mult_coef(n)
        end do

        call cpu_time(start)
        do l = 1, ms
            do n = 2 - mod(l, 2), ms, 2
                ix = mod(n + 1, 2)

                Int(n, l) = First(ix, n) * &
                        (Second(ix, l) * fcoef_middle(ix) + &
                                Second(ix + 2, l) * fcoef_upper(ix))

                do r = ix + 2, md - 2 * bucket - 1, 2 * bucket
                    do k = r, r + 2 * (bucket - 1), 2
                        value = First(k, n) * (Second(k - 2, l) * fcoef_lower(k) + &
                                Second(k, l) * fcoef_middle(k) + Second(k + 2, l) * fcoef_upper(k))
                        !                        write(*,*) 'value = ', value, 'first = ', first(k, n), 'second = ', second(k, l),&
                        !                                'mult_coef = ', mult_coef(k), 'func_coef = ', func_coef(m, k)
                        Int(n, l) = Int(n, l) + value
                        ! write(*,*) 'to ', n, ',', l, ' adding ', value
                    end do
                    if ((abs(value) < accuracy * abs(Int(n, l))).or.(abs(value) < border * abs(Int(1,1)))) then
                        !                         write(*,*) 'n = ', n, 'l = ', l, 'required r = ', r
                        exit
                    end if
                enddo
            enddo
        enddo
        call cpu_time(finish)
        ! write(*,*) 'time = ', finish - start
        !        call log_matrix(FILE_DESCRIPTOR(WARNING), 'Delta', Delta, .false., matrix_size)
        !write(*,*)

    end subroutine triple_dep_integral
end module integrals