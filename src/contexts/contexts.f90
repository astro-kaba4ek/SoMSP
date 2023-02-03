! Created by odawing on 17.08.22.

module contexts
    use regime
    use spheroidal_scatterer
    use wavelength_point
    use spheroidal
    use constants
    use regime
    use integrals
    use matrix
    use legendre_functions
    use utils
    use scattering_directions
    use angle
    use scattering_context_module
    use spheroidal_context_module
    use far_context_module
    use barber_context_module

    implicit none
    

    ! Global computation context contains 
    ! - the pointer to the scattering context
    ! - 1 item of every type of subcontext for mode transitions that require them
    ! each subcontext has a state, 0 is uncalculated state
    ! call to initialize resets all states of subcontexts
    ! 
    type :: ComputationContext
        type(ScatteringContext), pointer :: scattering_context
        type(SpheroidalContext), private :: spheroidal_context
        type(FarContext), private :: far_context
        type(BarberContext), private :: barber_context
    contains
        procedure :: initialize => initialize_computation_context
        ! procedure :: get_tmatrix_size_by_mode
        procedure :: get_spheroidal_context, get_far_context, get_barber_context
        procedure, private :: check_validity
    end type ComputationContext
contains

    subroutine check_validity(this)
        class(ComputationContext), intent(in) :: this
        
        if (.not. associated(this%scattering_context)) then
            write(LOG_FD, *) '{ERROR} scattering context was not initialized'
            call exit(1)
        endif
    end subroutine check_validity

    function get_spheroidal_context(this, target_state, m, fnum, lnum) result(res)
        class(ComputationContext), target, intent(inout) :: this
        integer, intent(in) :: target_state
        integer, intent(in) :: m, fnum, lnum

        type(SpheroidalContext), pointer :: res

        call this%check_validity()

        call this%spheroidal_context%initialize(m, fnum, lnum, this%scattering_context, target_state) 

        res => this%spheroidal_context

    end function get_spheroidal_context

    function get_far_context(this, m, spheroidal_lnum, spherical_lnum) result(res)
        class(ComputationContext), target, intent(inout) :: this
        integer, intent(in) :: m, spheroidal_lnum, spherical_lnum

        type(FarContext), pointer :: res

        call this%check_validity()

        call assert(this%spheroidal_context%state >= 1, '{ERROR} attempt to create far context before spheroidal')

        call assert(this%spheroidal_context%fnum >= max(spheroidal_lnum, spherical_lnum), &
                    '{ERROR} spheroidal context size is too small for far context') 

        call this%far_context%initialize(m, spheroidal_lnum, spherical_lnum, this%spheroidal_context%layers(0,1))

        res => this%far_context

    end function get_far_context

    function get_barber_context(this, m, lnum) result(res)
        class(ComputationContext), target, intent(inout) :: this
        integer, intent(in) :: m, lnum

        type(BarberContext), pointer :: res

        call this%check_validity()

        call this%barber_context%initialize(m, lnum, this%scattering_context%calculation_point%k)

        res => this%barber_context

    end function get_barber_context

    subroutine initialize_computation_context(this, m, lnum, spherical_lnum, scattering_context)
        class(ComputationContext), intent(out) :: this
        integer, intent(in) :: m, lnum, spherical_lnum
        type(ScatteringContext), target, intent(in) :: scattering_context

        this%scattering_context => scattering_context

        call this%spheroidal_context%reset()
        call this%far_context%reset()
        call this%barber_context%reset()

        call this%spheroidal_context%initialize(m, get_full_function_size(lnum, spherical_lnum), lnum, scattering_context, 1)

    end subroutine initialize_computation_context

end module contexts