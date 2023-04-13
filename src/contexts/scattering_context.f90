module scattering_context_module
    use regime
    use spheroidal_scatterer
    use wavelength_point
    use scattering_directions

    implicit none

    type :: ScatteringContext
        type(SpheroidalScatterer) :: scatterer
        type(ScatteringDirections) :: directions
        type(WavelengthPoint) :: calculation_point
    contains
        procedure :: initialize => initialize_scattering_context
    end type ScatteringContext

contains
    subroutine initialize_scattering_context(this, &
        f, nol, xv, ab, &
        alpha, &
        lambda, ri, &
        lnum, spherical_lnum, minm, maxm, &
        ntheta, theta0, theta1, nphi, phi0, phi1)
    class(ScatteringContext), intent(out) :: this
    integer, intent(in) :: f, nol, lnum, spherical_lnum, minm, maxm, ntheta, nphi
    real(knd), intent(in) :: xv(nol), ab(nol), alpha, lambda, theta0, theta1, phi0, phi1
    complex(knd) :: ri(0:nol)

    call this%scatterer%set(f, xv, ab, alpha, nol, lambda)
    call this%calculation_point%initialize(lambda, nol, ri)
    call this%directions%initialize(alpha, ntheta, theta0, theta1, nphi, phi0, phi1, minm, maxm)

    end subroutine initialize_scattering_context

end module scattering_context_module