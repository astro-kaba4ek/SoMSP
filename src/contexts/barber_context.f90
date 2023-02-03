module barber_context_module
    use regime
    use constants
    use legendre_functions
    use logging

    implicit none 

    type :: BarberContext
        integer :: m, lnum
        real(knd) :: k

        ! for spherical barber and mishch
        complex(knd), allocatable :: f(:,:), g(:), grev(:), xd(:), xdrev(:)
        ! for spherical mishch
        complex(knd), allocatable :: mishch_mult(:)

        integer :: state
    contains
        procedure, private :: check_barber
        procedure :: initialize => initialize_barber_context
        procedure :: reset => reset_barber_context
        final :: delete_barber_context
    end type BarberContext
contains
    subroutine reset_barber_context(this)
        class(BarberContext), intent(out) :: this

        this%state = 0
    end subroutine reset_barber_context

    subroutine initialize_barber_context(this, m, lnum, k)
        class(BarberContext), intent(inout) :: this
        integer, intent(in) :: m, lnum
        real(knd), intent(in) :: k

        if (.not. allocated(this%g) .or. m /= this%m .or. lnum > this%lnum) then
            this%state = 0
        endif

        call this%check_barber(m, lnum, k)
    end subroutine initialize_barber_context

    subroutine check_barber(this, m, lnum, k)
        class(BarberContext), intent(inout) :: this
        integer, intent(in) :: m, lnum
        real(knd), intent(in) :: k

        integer :: i
        real(knd) :: legendre_coef(lnum + 1)
        real(knd) :: start, finish

        if (this%state == 1) return

        this%m = m
        this%lnum = lnum
        this%k = k

        if (allocated(this%g) .and. (size(this%g) /=  lnum)) then
            deallocate(this%f, this%g, this%grev, this%xd, this%xdrev, this%mishch_mult)
        end if
        if (.not. allocated(this%g)) then
            allocate(this%f(lnum, lnum), &
                this%g(lnum), this%grev(lnum), &
                this%xd(lnum), this%xdrev(lnum), this%mishch_mult(lnum))
        end if

        this%f = 0
        this%g = 0
        this%grev = 0
        this%xd = 0
        this%xdrev = 0
        this%mishch_mult = 0

        legendre_coef = calculate_legendre_coef(m, lnum + 1)

        if (m == 0) then
            do i = 0, lnum - 1
                if (i > 0) then
                    this%f(i + 1, i) = 1.0_knd / sqrt((2.0_knd * i + 1) * (2.0_knd * i - 1))
                end if
                if (i < lnum - 1) then
                    ! m = 0
                    this%f(i + 1, i + 2) = 1.0_knd / sqrt((2.0_knd * i + 1) * (2.0_knd * i + 3))
                end if
            end do
        else
            do i = m, m + lnum - 1
                if (i > m) then
                    ! m > 0
                    this%f(i - m + 1, i - m) = sqrt(((i + m) * (i - m)) / &
                            (i ** 2 * (2.0_knd * i + 1) * (2.0_knd * i - 1)))
                end if
                if (i - m + 1 < lnum) then
                    ! m > 0
                    this%f(i - m + 1, i - m + 2) = sqrt(((i - m + 1.0_knd) * (i + m + 1.0_knd)) / &
                            ((2.0_knd * i + 1) * (2.0_knd * i + 3) * (i + 1) ** 2))
                end if
                ! 1/g1
                this%grev(i - m + 1) = i * (i + 1.0_knd) / m
                this%g(i - m + 1) = 1.0_knd / this%grev(i - m + 1)
            end do
        end if

        if (m == 0) then
            do i = 1, lnum - 1
                ! for m = 0
                this%xd(i + 1) = sqrt(2.0_knd / (2.0_knd * i + 1)) * (2.0_knd * i + 1) / (4 * i * (i + 1))
                this%xdrev(i + 1) = 1.0_knd / this%xd(i + 1)
            end do
        else
            do i = m, m + lnum - 1
                ! for m = 1
                !                xd(i, i) = sqrt(2q0 / (2q0 * i + 1) * (i + 1) * i) * (2q0 * i + 1) / (2 * (i * (i + 1))**2)
                this%xd(i - m + 1) = sqrt(legendre_coef(i - m + 1)) / (i * (i + 1.0_knd))
                !                xdrev(i, i) = 1q0 / (sqrt(2q0 / (2q0 * i + 1) * (i + 1) * i) * (2q0 * i + 1) / (2 * (i * (i + 1))**2))
                this%xdrev(i - m + 1) = 1.0_knd / this%xd(i - m + 1)
            end do
        end if

        if (this%m > 0) then
            this%mishch_mult = (/ (sqrt(real(i, knd) * (i + 1) / legendre_coef(i - m + 1)), i = m, m + lnum - 1) /)
        else
            ! this%mishch_mult(1) = 1
            this%mishch_mult = (/ (sqrt(real(i, knd) * (i + 1) / legendre_coef(i - m + 1)), i = 1, lnum) /)
        endif

        this%state = 1

        call cpu_time(finish)

        call log_time('context barber', finish - start)

    end subroutine check_barber

    subroutine delete_barber_context(this)
        type(BarberContext), intent(inout) :: this

        if (allocated(this%g)) then
            deallocate(this%f, this%g, this%grev, this%xd, this%xdrev, this%mishch_mult)
        end if

    end subroutine delete_barber_context

end module barber_context_module