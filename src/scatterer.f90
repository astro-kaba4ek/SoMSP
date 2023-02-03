! Created by drakosha on 15.07.2021.

module spheroidal_scatterer
    use regime
    use angle
    use constants
    use utils
    implicit none

    ! represents geometric properties of the scatterer
    type, public :: SpheroidalScatterer

        integer :: spheroidal_type

        !  parameters used during calculation obtained from shape

        integer :: number_of_layers
        !  radial coordinates of the scatterer's layers
        real(knd), allocatable, dimension(:) :: ksi   !  1..number_of_layers
        real(knd), allocatable, dimension(:) :: d   !  1..number_of_layers
        !  precalculated factor from extinction and scattering coefficient
        complex(knd), allocatable, dimension(:) :: c0   !  1..number_of_layers
        !  argument of the spheroidal function corresponding to the external field (layer 0)
        real(knd), allocatable, dimension(:) :: common_factor   !  1..number_of_layers

    contains

        procedure, public :: set

        final :: delete_scatterer
    end type SpheroidalScatterer

    type, public :: SpheroidalShape
        ! shape of the scatterrer

        integer :: spheroidal_type
        !  radius of the equivalent sphere
        real(knd) :: rv
        !  major semiaxis
        real(knd) :: a
        !  minor semiaxes
        real(knd) :: b
        !  a / b
        real(knd) :: ab

        type(AngleType) :: alpha
    contains
        procedure, public :: set => set_spheroidal_shape
    end type SpheroidalShape

contains
    !  for the given rv and ab returns a - the major semiaxis
    real(knd) function get_a(f, rv, ab) result(a)
        integer :: f
        real(knd) :: rv, ab

        if (f == 1) then
            a = rv * (ab ** (2q0 / 3q0))
        else
            a = rv * (ab ** (1q0 / 3q0))
        endif
    end function get_a

    !  for the given rv and ab returns b - the minor semiaxis
    real(knd) function get_b(f, rv, ab) result(b)
        integer :: f
        real(knd) :: rv, ab

        if (f == 1) then
            b = rv / ab ** (1q0 / 3q0)
        else
            b = rv / ab ** (2q0 / 3q0)
        endif
    end function get_b

    subroutine set(this, f, xv, ab, alpha, number_of_layers)

        class(SpheroidalScatterer) :: this
        real(knd) :: xv(number_of_layers), ab(number_of_layers), alpha, a(number_of_layers), b(number_of_layers)
        integer :: f, number_of_layers

        this%spheroidal_type = f

        if (allocated(this%ksi) .and. number_of_layers /= this%number_of_layers) then
            deallocate(this%ksi, this%c0, this%common_factor, this%d)
        end if
        if (.not. allocated(this%ksi)) then
            allocate(this%ksi(1:number_of_layers), this%c0(1:number_of_layers), this%common_factor(1:number_of_layers),&
            this%d(1:number_of_layers))
        end if
        this%number_of_layers = number_of_layers

        if (f == 1) then
            this%ksi = ab / sqrt(ab * ab - 1q0)
            this%c0 = (1q0 / ab)**(1q0 / 3q0)
        else
            this%ksi = 1q0 / sqrt(ab * ab - 1q0)
            this%c0 = (1q0 / ab)**(2q0 / 3q0)
        endif
        this%c0 = xv * sqrt(ab**2q0 - 1q0) * this%c0
        b = xv / ab**(1q0/3q0)
        a = ab * b
        this%d = (a * a - b * b) ** 0.5q0


        this%common_factor = 1q0 / (abs(this%c0)**2q0 * &
                sqrt((this%ksi**2 - this%spheroidal_type) * (this%ksi**2 - this%spheroidal_type * cos(alpha)**2)))

    end subroutine set

    subroutine set_spheroidal_shape(this, f, rv, ab, alpha)

        class(SpheroidalShape) :: this
        real(knd) :: rv, ab, alpha
        integer :: f

        this%spheroidal_type = f
        this%rv = rv
        this%a = get_a(f, rv, ab)
        this%b = this%a / ab
        this%ab = ab
        call this%alpha%set(alpha)

    end subroutine set_spheroidal_shape

    subroutine delete_scatterer(this)

        type(SpheroidalScatterer), intent(inout) :: this

        if (allocated(this%ksi)) then
            deallocate(this%ksi, this%c0, this%common_factor, this%d)
        end if

    end subroutine delete_scatterer

    real(knd) function convertQtoC(Q, f, rv, ab, alpha) result(C)
        integer, intent(in) :: f
        real(knd), intent(in) :: Q, rv, ab, alpha

        real(knd) :: a, b

        a = get_a(f, rv, ab)
        b = get_b(f, rv, ab)

        if (f == 1) then
            C = Q * PI * b * sqrt(a**2 * sin(alpha) ** 2 + b ** 2 * cos(alpha) ** 2)
        else
            C = Q * PI * a * sqrt(b**2 * sin(alpha) ** 2 + a ** 2 * cos(alpha) ** 2)
        end if
    end function convertQtoC

    type(ModeFactors) function get_c_factors_from_q(factors, shape)
        type(ModeFactors), intent(in) :: factors
        type(SpheroidalShape), intent(in) :: shape

        get_c_factors_from_q%Qext = convertQtoC(factors%Qext, shape%spheroidal_type, shape%rv, &
                shape%ab, shape%alpha%value)
        get_c_factors_from_q%Qsca = convertQtoC(factors%Qsca, shape%spheroidal_type, shape%rv, &
                shape%ab, shape%alpha%value)

    end function get_c_factors_from_q

    type(ModeFactors) function get_normalized_c_factors_from_q(factors, shape) result(result)
        type(ModeFactors), intent(in) :: factors
        type(SpheroidalShape), intent(in) :: shape

        result = get_c_factors_from_q(factors, shape)
        result%Qext = result%Qext / (PI * shape%rv ** 2)
        result%Qsca = result%Qsca / (PI * shape%rv ** 2)

    end function get_normalized_c_factors_from_q

end module spheroidal_scatterer