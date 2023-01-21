module scattering_directions
    use regime
    use constants
    use angle
    implicit none

    type :: ScatteringDirections
        !  angle between the spheroidal symmetry axis and the direction of incident light propagation
        type(AngleType) :: alpha
        ! scattering indicatrics is calculated for a square region [theta0..theta1]x[phi0..phi1] of (ntheta+1)*(nphi+1) points
        ! is skipped if ntheta*nphi = 0
        integer :: ntheta
        integer :: nphi

        type(AngleType), allocatable :: thetas(:)
        type(AngleType), allocatable :: phis(:)

    contains
        procedure :: initialize => initialize_scattering_directions
        procedure :: need_indicatrix
        final :: delete_scattering_directions
    end type ScatteringDirections

contains
    logical function need_indicatrix(this)
        class(ScatteringDirections), intent(in) :: this

        need_indicatrix = this%ntheta > 0 .and. this%nphi > 0
    end function need_indicatrix

    subroutine initialize_scattering_directions(this, alpha, ntheta, theta0, theta1, nphi, phi0, phi1, minm, maxm)
        class(ScatteringDirections), intent(out) :: this
        real(knd), intent(in) :: alpha, theta0, theta1, phi0, phi1
        integer, intent(in) :: ntheta, nphi, minm, maxm

        real(knd) :: dtheta, dphi, value
        integer :: i, j, k, primary_size, secondary_size, dt

        call this%alpha%set(alpha)
        this%ntheta = ntheta
        this%nphi = nphi

        if (.not. this%need_indicatrix()) then
            return
        endif

        if (nphi > 1) then
            dphi = (phi1 - phi0) / (nphi - 1)
        else
            dphi = 0
        endif

        this%phis = [(AngleType(phi0 + dphi * i), i = 0, nphi - 1)]

        if (ntheta > 1) then
            dtheta = (theta1 - theta0) / (ntheta - 1)
        else
            dtheta = 0
        endif

        if (allocated(this%thetas) .and. size(this%thetas) /= ntheta) deallocate(this%thetas)
        if (.not. allocated(this%thetas)) allocate(this%thetas(ntheta))

        do i = 0, ntheta - 1
            value = theta0 +  i * dtheta
            if (abs(value) < THETA_ADD .or. abs(PI - value) < THETA_ADD .or. abs(2 * PI - value) < THETA_ADD) then
                value = value + THETA_ADD
            endif
            this%thetas(i + 1) = AngleType(value)
        enddo

    end subroutine initialize_scattering_directions

    subroutine delete_scattering_directions(this)
        type(ScatteringDirections), intent(inout) :: this

        if (allocated(this%phis)) then
            deallocate(this%phis)
        endif

        if (allocated(this%thetas)) then
            deallocate(this%thetas)
        endif

    end subroutine delete_scattering_directions
end module scattering_directions