module mode_functors
    use regime
    use constants
    use contexts
    use spheroidal_initial
    use utils
    use spheroidal_scattering
    use tmatrix_conversion

implicit none


    abstract interface
        subroutine calculate_initial(computation_context, mode_item, initial)
            import :: knd
            import :: ComputationContext, ModeItem
            type(ComputationContext), intent(inout) :: computation_context
            type(ModeItem), intent(in) :: mode_item
            complex(knd), allocatable, intent(out) :: initial(:)
        end subroutine calculate_initial

        real(knd) function calculate_factor(computation_context, mode_item, solution)
            import :: knd
            import :: ComputationContext, ModeItem
            type(ComputationContext), intent(inout) :: computation_context
            type(ModeItem), intent(in) :: mode_item
            complex(knd), intent(in) :: solution(:)
        end function calculate_factor

    end interface

    type :: ModeFunctors
        procedure(calculate_initial), pointer, nopass :: initial_calculator
        procedure(calculate_factor), pointer, nopass :: extinction_calculator, scattering_calculator
    contains
    end type ModeFunctors

contains
    function get_mode_functions(mode) result(funcs)
        type(ModeInfo), intent(in) :: mode
        type(ModeFunctors) :: funcs

        if (mode == MODE_SPH_TE_UV) then
            funcs%initial_calculator => set_initial_uv_te
            funcs%extinction_calculator => get_extinction_factor_sph_uv
            funcs%scattering_calculator => get_scattering_factor_sph_uv
        elseif (mode == MODE_SPH_TM_UV) then
            funcs%initial_calculator => set_initial_uv_tm
            funcs%extinction_calculator => get_extinction_factor_sph_uv
            funcs%scattering_calculator => get_scattering_factor_sph_uv
        elseif (mode == MODE_SPH_TE_PQ) then
            funcs%initial_calculator => set_initial_pq_te
            funcs%extinction_calculator => get_extinction_factor_sph_pq
            funcs%scattering_calculator => get_scattering_factor_sph_pq
        elseif (mode == MODE_SPH_TM_PQ) then
            funcs%initial_calculator => set_initial_pq_tm
            funcs%extinction_calculator => get_extinction_factor_sph_pq
            funcs%scattering_calculator => get_scattering_factor_sph_pq
        elseif (mode == MODE_FAR_TE_UV) then
            funcs%initial_calculator => set_initial_far_uv_te
            funcs%extinction_calculator => get_extinction_factor_far_uv
            funcs%scattering_calculator => get_scattering_factor_far_uv
        elseif (mode == MODE_FAR_TM_UV) then
            funcs%initial_calculator => set_initial_far_uv_tm
            funcs%extinction_calculator => get_extinction_factor_far_uv
            funcs%scattering_calculator => get_scattering_factor_far_uv
        elseif (mode == MODE_FAR_TE_PQ) then
            funcs%initial_calculator => set_initial_far_pq_te
            funcs%extinction_calculator => get_extinction_factor_far_pq
            funcs%scattering_calculator => get_scattering_factor_far_pq
        elseif (mode == MODE_FAR_TM_PQ) then
            funcs%initial_calculator => set_initial_far_pq_tm
            funcs%extinction_calculator => get_extinction_factor_far_pq
            funcs%scattering_calculator => get_scattering_factor_far_pq
        else
            call assert(.false., 'can only calculate factors for spheroidal and far modes')
        endif
    end function get_mode_functions

end module mode_functors