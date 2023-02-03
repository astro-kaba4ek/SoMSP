module wavelength_point

    use regime
    use constants

    implicit none
    private

    ! represents a point of eps(lambda), mu(lambda) for calculation.
    type, public :: WavelengthPoint
        !  k = 2 * pi / lambda
        real(knd) :: k
        !  arrays from 0 to number of layers of the scatterer,
        !  where eps(0), mu(0) belong to the field outside of the scatterer
        complex(knd), allocatable, dimension(:) :: eps, mu

    contains
        procedure, private :: check_and_reallocate_arrays
        procedure, private :: initialize_from_eps_and_mu, initialize_from_refractive_index
        generic, public :: initialize => initialize_from_eps_and_mu, initialize_from_refractive_index
        procedure, public :: get_refractive_index
        final :: delete_wavelength_point
    end type WavelengthPoint

contains
    !  only (re)allocates eps and mu arrays if necessary
    subroutine check_and_reallocate_arrays(this, number_of_layers)
        class(WavelengthPoint), intent(inout) :: this
        integer, intent(in) :: number_of_layers

        if (allocated(this%eps) .and. size(this%eps) /= number_of_layers + 1) then
            deallocate(this%eps, this%mu)
        end if

        if (.not.allocated(this%eps)) then
            allocate(this%eps(0:number_of_layers), this%mu(0:number_of_layers))
        end if

    end subroutine check_and_reallocate_arrays

    complex(knd) function get_refractive_index(this, layer)
        class(WavelengthPoint), intent(in) :: this
        integer, intent(in) :: layer

        get_refractive_index = sqrt(this%eps(layer) * this%mu(layer))
    end function get_refractive_index

    !  base constructor from arrays [0:number_of_layers], 0 is the outside medium
    !  logs the point to info
    subroutine initialize_from_eps_and_mu(this, lambda, number_of_layers, eps, mu)

        class(WavelengthPoint), intent(inout) :: this
        real(knd), intent(in) :: lambda
        integer, intent(in) :: number_of_layers
        complex(knd), intent(in) :: eps(0:number_of_layers), mu(0:number_of_layers)

        this%k = 2q0 * PI / lambda

        call this%check_and_reallocate_arrays(number_of_layers)

        this%eps = eps
        this%mu = mu

    end subroutine initialize_from_eps_and_mu

    !  constructor from the array ri[0:number_of_layers]. assumes
    !  eps = ri**2
    !  mu = 1
    subroutine initialize_from_refractive_index(this, lambda, number_of_layers, ri)

        class(WavelengthPoint), intent(inout) :: this
        real(knd), intent(in) :: lambda
        integer, intent(in) :: number_of_layers
        complex(knd), intent(in) :: ri(0:number_of_layers)

        integer :: i

        call this%initialize_from_eps_and_mu(lambda, number_of_layers, ri**2, &
                (/(cmplx(1q0, 0q0, knd), i=0,number_of_layers)/))

    end subroutine initialize_from_refractive_index

    subroutine delete_wavelength_point(this)
        type(WavelengthPoint), intent(inout) :: this

        if (allocated(this%eps)) then
            deallocate(this%eps, this%mu)
        end if
    end subroutine delete_wavelength_point

end module wavelength_point