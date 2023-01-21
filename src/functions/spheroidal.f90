!  Wrapper for subroutines for spheroidal function calculation by Anri van Buren
!  from modules vb_prolate and vb_oblate
!  the kind parameter for complex function values is set in module regime as knd
module spheroidal

    use regime
    use constants

    use complex_prolate_swf
    use complex_oblate_swf
    use param
    use logging

    implicit none
    private
    public :: get_ksi, get_c, theoretical_wronskian_radial

    real(knd), parameter :: lower_bound = 1.0e-128_knd
    integer, parameter :: lower_exp = -300

    !  Represents the spheroidal function
    !  contains the tabulated radial functions of 1st, 2nd, 3rd kind and their derivatives, R^(i)_{m, n}(c, \ksi)
    !  angular functions of the 1st kind and it's derivative, S^(1)_{m, n}(c, \eta)
    !  and legendre coefficients for the angular function
    !  normalized with Mexxie-Shafke narmalization scheme to 1
    !  its calculate method calls subroutines calculating spheroidal functions by Anri van Buren:
    !  cprofcn (module vb_prolate) and coblfcn (module vb_oblate)
    type, public :: SpheroidalCalculation
        !  spheroidal_type = 1 for prolate s.f. and -1 for oblate s.f.
        integer :: spheroidal_type
        ! spheroidal function parameter m
        integer :: m
        !  lnum = n - m + 1 >= 0
        integer :: lnum
        ! spheroidal function parameter c
        complex(knd) :: c

        !  ksi = \ksi - 1 for prolate functions and
        !  ksi = \ksi for oblate functions
        !  so ksi \in [0..)
        real(knd) :: ksi
        ! radial functions and their derivatives for n - m + 1 = 1,.., lnum, arrays [1:lnum]
        complex(knd), allocatable, dimension(:) :: r1, r1d, r2, r2d, r3, r3d

        !  number of values of \eta for angular function
        integer :: narg
        !  array[1:narg] with parameters \eta of angular function
        real(knd), allocatable, dimension(:) :: arg
        !  angular function and its derivative for \eta \in arg, n - m = 0,.., lnum - 1
        !  arrays[1:lnum][1:narg]
        complex(knd), allocatable, dimension(:, :) :: s1, s1d

        !  number of legendre coefficients calculated (including zeros)
        integer :: maxd
        !  legendre coefficients for angular functions, arrays [0:maxd][1:lnum]
        complex(knd), allocatable, dimension(:, :) :: legendre

        logical :: functions_calculated
        logical :: legendre_coefficients_calculated
    contains
        !  Preparation
        !  (re)allocates arrays of function values to required sizes if necessary
        procedure, private :: check_and_reallocate_arrays

        !  Calculation
        !  calculates spheroidal functions with giver parameters
        procedure, public :: calculate => calculate_spheroidal_functions

        !  get array of wronskian values
        procedure :: get_wronskian

        !  Free memory if it was allocated (noexception)
        procedure :: delete_radial_functions, delete_angular_functions, delete_legendre_coefficients, &
                delete_calculation_full

        !  destructor, calls delete_calculation_full
        final :: delete_calculation
    end type SpheroidalCalculation

contains
    subroutine set(this, m, n, c, ksi, narg, arg, f)
        class(SpheroidalCalculation) :: this
        integer, intent(in) :: m, n, narg
        real(knd), intent(in) :: ksi, arg(narg)
        complex(knd), intent(in) :: c
        integer, optional, intent(in) :: f

        if (present(f)) then
            this%spheroidal_type = f
        else
            this%spheroidal_type = 1
        endif

        this%m = m
        this%c = c
        if (this%spheroidal_type == 1) then
            this%ksi = ksi - 1q0
        else
            this%ksi = ksi
        endif
        this%narg = narg
        this%lnum = n - m + 1

        if (.not.(allocated(this%r1) .and. (n - m + 1 == this%lnum))) then
            !write(*,*) 'allocating spheroidal calculation'
            if (allocated(this%r1)) then
                call delete_calculation(this)
            endif

            this%lnum = n - m + 1

            allocate(this%r1(this%lnum), this%r1d(this%lnum), &
                    this%r2(this%lnum), this%r2d(this%lnum), this%r3(this%lnum), this%r3d(this%lnum))
            allocate(this%s1(this%lnum, this%narg), this%s1d(this%lnum, this%narg))
            allocate(this%arg(narg))
        endif

        if (.not.(allocated(this%arg).and. (size(this%arg) == narg))) then
            call this%delete_angular_functions()
            allocate(this%s1(this%lnum, this%narg), this%s1d(this%lnum, this%narg))
            allocate(this%arg(narg))
        endif
        !allocate(this%arg(this%narg))
        this%arg = arg
        this%functions_calculated = .false.

    end subroutine set


    ! checks if the arrays for function values and angular function argument values need to be reallocated
    ! assignes and reallocates them only if necessary
    subroutine check_and_reallocate_arrays(this, lnum, narg, allocate_radial, allocate_angular)
        class(SpheroidalCalculation), intent(inout) :: this
        integer, intent(in) :: lnum, narg
        logical, intent(in) :: allocate_radial, allocate_angular

        if (allocate_radial .and. allocated(this%r1) .and. lnum /= this%lnum) then
            deallocate(this%r1, this%r1d, this%r2, this%r2d, this%r3, this%r3d, this%s1, this%s1d)
        end if

        if (allocated(this%arg) .and. narg /= this%narg) then
            deallocate(this%arg)
        end if

        if (allocate_angular .and. allocated(this%s1) .and. (narg /= this%narg .or. lnum /= this%lnum)) then
            deallocate(this%s1, this%s1d)
        end if

        this%lnum = lnum
        this%narg = narg
        if (allocate_radial .and. .not. allocated(this%r1)) then
            allocate(this%r1(lnum), this%r1d(lnum), this%r2(lnum), this%r2d(lnum), this%r3(lnum), this%r3d(lnum))
        end if

        if (.not. allocated(this%arg)) then
            allocate(this%arg(narg))
        end if
        if (allocate_angular .and. .not. allocated(this%s1)) then
            allocate(this%s1(lnum, narg), this%s1d(lnum, narg))
        end if
    end subroutine check_and_reallocate_arrays

    subroutine calculate_spheroidal_functions(this, m, lnum, c, ksi, narg, arg, f, &
            need_radial_functions, need_angular_functions, need_legendre_coefficients)
        class(SpheroidalCalculation), intent(inout) :: this
        integer, intent(in) :: m, lnum, narg, f
        real(knd), intent(in) :: ksi, arg(narg)
        complex(knd), intent(in) :: c
        logical, optional, intent(in) :: need_radial_functions, need_angular_functions, need_legendre_coefficients
        integer :: legendre_exp(0:MAXIMUM_D_COEFF_NUMBER), maxb, acts, md
        complex(knd) :: value
        real(knd) :: rad, start, finish, middle

        !  cprofcn and coblfcn calculate exponent and mantissa in different arrays
        !  these arrays are used to store them, their types have kind = vb_kind
        !  the arrays from SpheroidalCalculation contain the full value
        complex(vb_knd), allocatable, dimension(:) :: r1, r1d, r2, r2d
        integer, allocatable, dimension(:) :: r1_exp, r1d_exp, r2_exp, r2d_exp
        complex(vb_knd), allocatable, dimension(:, :) :: s1, s1d
        integer, allocatable, dimension(:, :) :: s1_exp, s1d_exp
        !  legendre coeffictient ratios,
        !  the array is allocated inside cprofcn, coblfcn and deallocated right after
        !  this%legendre is initialized from it
        integer :: maxd, mmaxval
        complex(vb_knd), allocatable, dimension(:, :) :: enr
        !  expected accuracies for function values
        integer, allocatable, dimension(:) :: naccr
        integer, allocatable, dimension(:, :) :: naccs, naccds
        real(vb_knd) :: vb_arg(narg)

        integer :: i, j, mode_radial, mode_angular
        logical :: need_radial, need_angular, need_legendre

        complex(knd), allocatable :: small_leg(:,:)

!        complex(knd), allocatable :: legtmp(:,:)

        !  process optional parameters. used for indicatrix calculation
        !  where only angular function values are needed
        need_radial = .true.
        need_angular = .true.
        need_legendre = .true.
        if (present(need_radial_functions)) then
            need_radial = need_radial_functions
        end if
        if (present(need_angular_functions)) then
            need_angular = need_angular_functions
        end if
        if (present(need_legendre_coefficients)) then
            need_legendre = need_legendre_coefficients
        end if
        if (need_radial) then
            mode_radial = 2
        else
            mode_radial = 0
        end if
        if (need_angular) then
            mode_angular = 2
        else
            mode_angular = 0
        end if

        allocate(r1(lnum), r1_exp(lnum), r1d(lnum), r1d_exp(lnum), &
                r2(lnum), r2_exp(lnum), r2d(lnum), r2d_exp(lnum))
        r1 = 0; r1d = 0; r2 = 0; r2d = 0; r1_exp = 0; r1d_exp = 0; r2_exp = 0; r2d_exp = 0

        allocate(s1(lnum, narg), s1_exp(lnum, narg), s1d(lnum, narg), s1d_exp(lnum, narg))
        s1 = 0; s1d = 0; s1_exp = 0; s1d_exp = 0

        allocate(naccr(lnum), naccs(lnum, narg), naccds(lnum, narg))
        call cpu_time(start)
        !  call functions by van Buren
        !  use different arrays to support different kind parameters for van Buren's
        !  procedures and our procedures
        if (f == 1) then
            call cprofcn(cmplx(c%re, c%im, vb_knd), m, lnum, mode_radial, real(ksi, vb_knd) - 1, &
                    mode_angular, 1, narg, real(arg, vb_knd), &
                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, naccr, &
                    s1, s1_exp, s1d, s1d_exp, naccs, naccds, need_legendre, maxd, enr)
        else
            call coblfcn(cmplx(c%re, c%im, vb_knd), m, lnum, mode_radial, real(ksi, vb_knd), &
                    mode_angular, 1, narg, real(arg, vb_knd), &
                    r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp, naccr, &
                    s1, s1_exp, s1d, s1d_exp, naccs, naccds, need_legendre, maxd, enr)
        end if

        !write(*,*) 'enr1 = ', enr(1:1000, 1)
        !write(*,*) 's1_base = ', s1(1:10, 1)

        !  set arguments
        this%spheroidal_type = f
        this%m = m
        this%c = c
        this%ksi = ksi
        call this%check_and_reallocate_arrays(lnum, narg, need_radial, need_angular)
        this%arg = arg
        if (LOG_INFO) write(LOG_FD, *) '{INFO} calculated functions for m = ', this%m, ' lnum = ', this%lnum, ' c = ', this%c, &
            ' ksi = ', this%ksi, ' nu = ', this%arg

        !write(*,*) 'r1 = ', r1(1)
        !  set function values
        if (need_radial) then
            this%r1 = 0
            this%r1d = 0
            this%r2 = 0
            this%r2d = 0
            this%r3 = 0
            this%r3d = 0

            this%r1 = r1 * (10q0**r1_exp)
            this%r1d = r1d * (10q0**r1d_exp)
            this%r2 = r2 * (10q0**r2_exp)
            this%r2d = r2d * (10q0**r2d_exp)
            this%r3 = this%r1 + cmplx(0q0, 1q0, knd) * this%r2
            this%r3d = this%r1d + cmplx(0q0, 1q0, knd) * this%r2d
        end if
        if (need_angular) then
            this%s1 = 0
            this%s1d = 0
            this%s1 = s1 * (10q0**s1_exp)
            this%s1d = s1d * (10q0**s1d_exp)
        end if
        this%functions_calculated = .true.

        !call check_expected_accuracy(lnum, narg, naccr, naccs, naccds)
        deallocate(r1, r1_exp, r1d, r1d_exp, r2, r2_exp, r2d, r2d_exp)
        deallocate(s1, s1_exp, s1d, s1d_exp)
        deallocate(naccr, naccs, naccds)
        call cpu_time(finish)
        call log_time('function values', finish - start)
        ! write(*,*) 'finctions time = ', finish - start
        !  set legendre coefficients
        call cpu_time(start)
        if (need_legendre) then
            !  the enr array obtained from cprofcn and coblfcn has the ratios of
            !  legendre coefficients enr(i) = legendre(2 * i + 2 + ix) / legendre(2 * i + ix)
            !  ix = (n - m) mod 2
            if (allocated(this%legendre)) then
                deallocate(this%legendre)
            endif

            this%maxd = min(2 * maxd + 1, MAXIMUM_D_COEFF_NUMBER)
            allocate(this%legendre(0:this%maxd, this%lnum))
            this%legendre = 0
            rad = radix(real(enr(1,1)))
!            write(*,*) 'rad = ', rad
            mmaxval = 0
            md = 0
            do i = 1, this%lnum
                legendre_exp = 0
                this%legendre(mod(i + 1, 2), i) = 1.0q0
                !write(*,*) 'i = ', i, 'enr = ', enr(:100,i)
                acts = this%maxd
                do j = mod(abs(i + 1), 2) + 2, this%maxd, 2
                    value = this%legendre(j - 2, i) * enr(j / 2, i)
!                    write(*,*) 'real_exp = ', exponent(real(value, knd)), 'imag_exp = ', exponent(imag(value))
                    legendre_exp(j) = exponent(real(value, knd)) + legendre_exp(j - 2)
                    this%legendre(j, i) = value / rad ** exponent(real(value, knd))
                    mmaxval = max(mmaxval, legendre_exp(j))
                    if ((legendre_exp(j) < lower_exp) .or. (abs(value * rad ** legendre_exp(j)) < lower_bound)) then
                        acts = j
                        exit
                    end if
                enddo
                do j = mod(abs(i + 1), 2), acts, 2
                    legendre_exp(j) = legendre_exp(j) - mmaxval
                this%legendre(j,i) = this%legendre(j,i) * rad ** legendre_exp(j)

                    end do
                md = max(md, acts)
            enddo

            deallocate(enr)
            call cpu_time(middle)
            md = max(md, lnum * 2)
            this%maxd = md
            allocate(small_leg(0:md,lnum))
            small_leg(0:md,1:lnum) = this%legendre(0:md, 1:lnum)
            deallocate(this%legendre)
            allocate(this%legendre(0:md,1:lnum))
            this%legendre = small_leg
            deallocate(small_leg)
                call normalize(this%legendre, md, this%m)

            this%legendre_coefficients_calculated = .true.
        else
            if (allocated(enr)) then
                deallocate(enr)
            endif
        endif
        call cpu_time(finish)
        call log_time('legendre coefficients', finish - start)

    end subroutine calculate_spheroidal_functions

    function get_wronskian(this) result(WR)
        class(SpheroidalCalculation) :: this
        complex(knd) :: WR(this%lnum)

        !write(*, *) this%r1 * this%r2d
        !write(*, *) this%r2 * this%r1d
        !write(*,*) this%r1 * this%r2d -
        WR = this%r1 * this%r2d - this%r2 * this%r1d
    end function get_wronskian

    complex(knd) function theoretical_wronskian_radial(c, ksi, f) result(wr)
        complex(knd) :: c
        real(knd) :: ksi
        integer :: f
        if (f == 1) then
            wr = 1.0q0 / c / (ksi * ksi - 1.0q0)
        else
            wr = 1.0q0 / c / (ksi * ksi + 1.0q0)
        end if
    end function theoretical_wronskian_radial


    complex(knd) function get_c(xv, ab, f, ri)
        real(knd) :: xv, ab
        integer :: f
        complex(knd), optional :: ri

        if (.not. present(ri)) then
            ri = 1q0
        end if

        if (f == 1) then
            get_c = (1q0 / ab)**(1q0 / 3q0)
        else
            get_c = (1q0 / ab)**(2q0 / 3q0)
        endif
        get_c = xv * sqrt(ab**2q0 - 1q0) * get_c * ri

    end function get_c

    real(knd) function get_ksi(ab, f)
        real(knd) :: ab
        integer :: f

        if (f == 1) then
            get_ksi = ab / sqrt(ab * ab - 1q0)
        else
            get_ksi = 1q0 / sqrt(ab * ab - 1q0)
        endif

    end function get_ksi


    subroutine delete_calculation(this)
        type(SpheroidalCalculation), intent(inout) :: this

        call this%delete_calculation_full()

    end subroutine delete_calculation

    subroutine delete_radial_functions(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%functions_calculated = .false.
        if (allocated(this%r1)) then
            deallocate(this%r1, this%r1d, this%r2, this%r2d, this%r3, this%r3d)
        endif

    end subroutine delete_radial_functions

    subroutine delete_angular_functions(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%functions_calculated = .false.
        if (allocated(this%arg)) then
            deallocate(this%arg, this%s1, this%s1d)
        endif

    end subroutine delete_angular_functions

    subroutine delete_legendre_coefficients(this)
        class(SpheroidalCalculation), intent(inout) :: this

        this%legendre_coefficients_calculated = .false.
        if (allocated(this%legendre)) then
            deallocate(this%legendre)
        endif

    end subroutine delete_legendre_coefficients

    subroutine delete_calculation_full(this)
        class(SpheroidalCalculation), intent(inout) :: this

        call this%delete_radial_functions()
        call this%delete_angular_functions()
        call this%delete_legendre_coefficients()

    end subroutine delete_calculation_full

    !  normalizes the array of legendre coefficients d of the size maxd so that
    !  \sum_{r=0}^\infty d[r]^2 * 2/(2r+2m+1) * (r+2m)!/r! = 1
    subroutine normalize(d, maxd, m)
        integer :: n, r, maxd, m, j
        complex(knd) d(0:, :), norm, coef(0:maxd)
        real(knd) :: factorial

        factorial = 1q0
        do r = 2, 2 * m
            factorial = factorial * r
        end do
        do r = 0, maxd
            coef(r) = factorial* 2q0 / (2 * r + 2 * m + 1)
            factorial = factorial * (r + 2 * m + 1) / (r + 1)
        end do
        ! write(*,*) size(d(0,:)), size(d(:,1)), maxd
        do j = 1, size(d(0,:))
        norm = 0q0
        do r = 0, maxd
            norm = norm + d(r, j) ** 2 * coef(r)
        enddo
        !write(*,*) 'norm = ', norm
        norm = 1.0_knd / sqrt(norm)
        d(:maxd, j) = d(:maxd, j) * norm

        end do
    end subroutine normalize

end module spheroidal