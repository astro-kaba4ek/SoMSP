module calculation_models
    use regime
    use utils
    implicit none
    private

    type, public :: SolutionForIndicatrix
        integer :: m
        integer :: basis_type
        integer :: source_tm
        integer :: source_te
    end type SolutionForIndicatrix

    type, abstract, public :: CalculationModel
    contains
        procedure(build_queue), deferred, nopass :: build_mode_queue
        procedure(print_row), deferred, nopass :: print_mode_row
        procedure(get_solution_places), deferred, nopass :: solution_places
    end type CalculationModel

    abstract interface
        subroutine build_queue(m, lnum, spherical_lnum, queue)
            import Node
            type(Node), allocatable, intent(out) :: queue(:)
            integer, intent(in) :: m, lnum, spherical_lnum
        end subroutine build_queue


        subroutine print_row(m, mode_res)
            import ModeCalculationResult
            integer, intent(in) :: m
            type(ModeCalculationResult), intent(in) :: mode_res(:)
        end subroutine print_row

        function get_solution_places(m) result(res)
            import :: SolutionForIndicatrix

            integer, intent(in) :: m
            type(SolutionForIndicatrix), allocatable :: res(:)
        end function get_solution_places
    end interface

    type, extends(CalculationModel) :: AllUVModel
    contains
        procedure, nopass :: build_mode_queue => build_all_uv_queue
        procedure, nopass :: print_mode_row => print_all_uv_table_row
        procedure, nopass :: solution_places => solution_places_all_uv
    end type AllUVModel

    type, extends(CalculationModel) :: CompareTEModel
    contains
        procedure, nopass :: build_mode_queue => build_queue_compare_te
        procedure, nopass :: print_mode_row => print_row_compare_te
        procedure, nopass :: solution_places => solution_places_compare_te
    end type CompareTEModel

    type, extends(CalculationModel) :: CompareUVPQModel
    contains
        procedure, nopass :: build_mode_queue => build_queue_compare_uv_pq
        procedure, nopass :: print_mode_row => print_row_compare_uv_pq
        procedure, nopass :: solution_places => solution_places_compare_uv_pq
    end type CompareUVPQModel

    type, extends(CalculationModel) :: UVPQModel
    contains
        procedure, nopass :: build_mode_queue => build_queue_uv_pq
        procedure, nopass :: print_mode_row => print_row_uv_pq
        procedure, nopass :: solution_places => solution_places_uv_pq
    end type UVPQModel

    type, extends(CalculationModel) :: UVPQTEFromTMModel
    contains
        procedure, nopass :: build_mode_queue => build_queue_uv_pq_te_from_tm
        procedure, nopass :: print_mode_row => print_row_uv_pq_te_from_tm
        procedure, nopass :: solution_places => solution_places_uv_pq_te_from_tm
    end type UVPQTEFromTMModel

    type, extends(CalculationModel) :: UVPQTEFromTMWithFarModel
    contains
        procedure, nopass :: build_mode_queue => build_queue_uv_pq_te_from_tm_with_far
        procedure, nopass :: print_mode_row => print_row_uv_pq_te_from_tm_with_far
        procedure, nopass :: solution_places => solution_places_uv_pq_te_from_tm
    end type UVPQTEFromTMWithFarModel

    interface CalculationModel
        procedure :: construct_calculation_model
    end interface CalculationModel

    character(*), parameter :: row_format = '(1x,1I5,1x,1A12,1x,6E24.15)'
    ! 100 format(' ',1I5, ' ', 1A12,' ',6E24.15)
contains
    subroutine build_queue_uv_pq_te_from_tm(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        if (m == 1) then
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .false., .false.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TE_UV, m, spherical_lnum, [3], .false., .false.), &
                Node(MODE_SPH_TE_UV, m, lnum, [4], .true., .true.), &
                Node(MODE_SPH_TM_PQ, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_PQ, m, lnum, [integer::], .true., .true.) &
            ]
        elseif (m == 0) then
            queue = [Node::]
        else
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .false., .false.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TE_UV, m, spherical_lnum, [3], .false., .false.), &
                Node(MODE_SPH_TE_UV, m, lnum, [4], .true., .true.) &
            ]
        endif

    end subroutine build_queue_uv_pq_te_from_tm

    subroutine build_queue_uv_pq_te_from_tm_with_far(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        if (m == 1) then
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .true., .true.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TE_UV, m, spherical_lnum, [3], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [4], .true., .true.), &
                Node(MODE_SPH_TM_PQ, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_PQ, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_PQ, m, lnum, [6], .true., .true.), &
                Node(MODE_FAR_TE_PQ, m, lnum, [7], .true., .true.) &
            ]
        elseif (m == 0) then
            queue = [Node::]
        else
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .true., .true.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TE_UV, m, spherical_lnum, [3], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [4], .true., .true.) &
            ]
        endif

    end subroutine build_queue_uv_pq_te_from_tm_with_far

    subroutine print_row_uv_pq_te_from_tm(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        if (m > 0) then
            write(*,row_format) m, 'UV', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(5)%factors%Qext, &
                mode_res(5)%factors%Qsca, &
                mode_res(5)%factors%qabs()
        endif
        if (m == 1) then
            write(*,row_format) m, 'PQ', &
                mode_res(6)%factors%Qext, &
                mode_res(6)%factors%Qsca, &
                mode_res(6)%factors%qabs(), &
                mode_res(7)%factors%Qext, &
                mode_res(7)%factors%Qsca , &
                mode_res(7)%factors%qabs()  
        endif
    end subroutine print_row_uv_pq_te_from_tm

    subroutine print_row_uv_pq_te_from_tm_with_far(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        if (m > 0) then
            write(*,row_format) m, 'UV', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(5)%factors%Qext, &
                mode_res(5)%factors%Qsca, &
                mode_res(5)%factors%qabs()
            write(*,row_format) m, 'FAR_UV', &
                mode_res(2)%factors%Qext, &
                mode_res(2)%factors%Qsca, &
                mode_res(2)%factors%qabs(), &
                mode_res(4)%factors%Qext, &
                mode_res(4)%factors%Qsca, &
                mode_res(4)%factors%qabs()
        endif
        if (m == 1) then
            write(*,row_format) m, 'PQ', &
                mode_res(6)%factors%Qext, &
                mode_res(6)%factors%Qsca, &
                mode_res(6)%factors%qabs(), &
                mode_res(7)%factors%Qext, &
                mode_res(7)%factors%Qsca , &
                mode_res(7)%factors%qabs()  
            write(*,row_format) m, 'FAR_PQ', &
                mode_res(8)%factors%Qext, &
                mode_res(8)%factors%Qsca, &
                mode_res(8)%factors%qabs(), &
                mode_res(9)%factors%Qext, &
                mode_res(9)%factors%Qsca , &
                mode_res(9)%factors%qabs()  
        endif
    end subroutine print_row_uv_pq_te_from_tm_with_far

    function solution_places_uv_pq_te_from_tm(m) result(res)
        integer, intent(in) :: m
        type(SolutionForIndicatrix), allocatable :: res(:)

        if (m == 0) then
            res = [SolutionForIndicatrix::]
        elseif(m == 1) then
            res = [&
                SolutionForIndicatrix(1, PQ, 6, 7), &
                SolutionForIndicatrix(1, UV, 1, 5) &
            ]
        else
            res = [SolutionForIndicatrix(m, UV, 1, 5)]
        endif
    end function solution_places_uv_pq_te_from_tm

    subroutine build_queue_uv_pq(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        if (m == 1) then
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TM_PQ, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_PQ, m, lnum, [integer::], .true., .true.) &
            ]
        elseif (m == 0) then
            queue = [Node::]
        else
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [integer::], .true., .true.) &
            ]     
        endif

    end subroutine build_queue_uv_pq

    subroutine print_row_uv_pq(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        if (m > 0) then
            write(*,row_format) m, 'UV', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(2)%factors%Qext, &
                mode_res(2)%factors%Qsca, &
                mode_res(2)%factors%qabs()
        endif
        if (m == 1) then
            write(*,row_format) m, 'PQ', &
                mode_res(3)%factors%Qext, &
                mode_res(3)%factors%Qsca, &
                mode_res(3)%factors%qabs(), &
                mode_res(4)%factors%Qext, &
                mode_res(4)%factors%Qsca, &
                mode_res(4)%factors%qabs()   
        endif
    end subroutine print_row_uv_pq

    function solution_places_uv_pq(m) result(res)
        integer, intent(in) :: m
        type(SolutionForIndicatrix), allocatable :: res(:)

        if (m == 0) then
            res = [SolutionForIndicatrix::]
        elseif(m == 1) then
            res = [&
                SolutionForIndicatrix(1, PQ, 3, 4), &
                SolutionForIndicatrix(1, UV, 1, 2) &
            ]
        else
            res = [SolutionForIndicatrix(m, UV, 1, 2)]
        endif
    end function solution_places_uv_pq

    subroutine build_queue_compare_te(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        if (m == 0) then
            allocate(queue(0))
        else
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .false., .false.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TE_UV, m, spherical_lnum, [3], .false., .false.), &
                Node(MODE_SPH_TE_UV, m, spherical_lnum, [4], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [integer::], .true., .false.) &
            ]
        endif

    end subroutine build_queue_compare_te

    subroutine print_row_compare_te(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        if (m > 0) then
            write(*,row_format) m, 'UV', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(6)%factors%Qext, &
                mode_res(6)%factors%Qsca, &
                mode_res(6)%factors%qabs()
            write(*,row_format) m, 'TE_from_TM', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(5)%factors%Qext, &
                mode_res(5)%factors%Qsca, &
                mode_res(5)%factors%qabs()
        endif
    end subroutine print_row_compare_te

    function solution_places_compare_te(m) result(res)
        integer, intent(in) :: m
        type(SolutionForIndicatrix), allocatable :: res(:)

        if (m == 0) then
            res = [SolutionForIndicatrix::]
        else
            res = [SolutionForIndicatrix(m, UV, 1, 5)]
        endif
    end function solution_places_compare_te

    subroutine build_queue_compare_uv_pq(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        if (m == 0) then
            queue = [ &
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_UV, m, lnum, [integer::], .true., .true.) &
            ]
        elseif (m == 1) then
            queue = [ &
                Node(MODE_SPH_TM_PQ, m, lnum, [integer::], .true., .true.), &
                Node(MODE_SPH_TE_PQ, m, lnum, [integer::], .true., .true.) &
            ]
        else 
            queue = [Node::]
        endif

    end subroutine build_queue_compare_uv_pq

    subroutine print_row_compare_uv_pq(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        if (m == 0) then
            write(*,row_format) m, 'UV', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(2)%factors%Qext, &
                mode_res(2)%factors%Qsca, &
                mode_res(2)%factors%qabs()
        elseif (m == 1) then
            write(*,row_format) m, 'PQ', &
                mode_res(1)%factors%Qext, &
                mode_res(1)%factors%Qsca, &
                mode_res(1)%factors%qabs(), &
                mode_res(2)%factors%Qext, &
                mode_res(2)%factors%Qsca, &
                mode_res(2)%factors%qabs()
        endif
    end subroutine print_row_compare_uv_pq

    function solution_places_compare_uv_pq(m) result(res)
        integer, intent(in) :: m
        type(SolutionForIndicatrix), allocatable :: res(:)

        if (m == 0) then
            res = [SolutionForIndicatrix(m, PQ, 1, 2)]
        elseif (m == 1) then
            res = [SolutionForIndicatrix(m, UV, 1, 2)]
        else
            res = [SolutionForIndicatrix::]
        endif
    end function solution_places_compare_uv_pq

    integer function unimplemented_func(m)
        integer, intent(in) :: m

        call assert(.false., 'function is not implemented')
    end function unimplemented_func
    
    subroutine build_all_uv_queue(m, lnum, spherical_lnum, queue)
        type(Node), allocatable, intent(out) :: queue(:)
        integer, intent(in) :: m, lnum, spherical_lnum

        queue = [ &
            Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .true.), &
            Node(MODE_SPH_TE_UV, m, lnum, [integer::], .true., .true.) &
        ]

    end subroutine build_all_uv_queue

    subroutine print_all_uv_table_row(m, mode_res)
        integer, intent(in) :: m
        type(ModeCalculationResult), intent(in) :: mode_res(:)

        write(*,row_format) m, 'UV', &
            mode_res(1)%factors%Qext, &
            mode_res(1)%factors%Qsca, &
            mode_res(1)%factors%qabs(), &
            mode_res(2)%factors%Qext, &
            mode_res(2)%factors%Qsca, &
            mode_res(2)%factors%qabs()
    end subroutine print_all_uv_table_row

    function solution_places_all_uv(m) result(res)
        integer, intent(in) :: m
        type(SolutionForIndicatrix), allocatable :: res(:)

        res = [SolutionForIndicatrix(m, UV, 1, 2)]

    end function solution_places_all_uv


    function construct_calculation_model(model_name) result(model)
        character(*), intent(in) :: model_name
        class(CalculationModel), allocatable :: model

        if (model_name == 'all_uv') then
            allocate(AllUVModel::model)
        elseif (model_name == 'compare_te') then
            allocate(CompareTEModel::model)
        elseif (model_name == 'compare_uv_pq') then
            allocate(CompareUVPQModel::model)
        elseif (model_name == 'uv_pq') then
            allocate(UVPQModel::model)
        elseif (model_name == 'uv_pq_te_from_tm') then
            allocate(UVPQTEFromTMModel::model)
        elseif (model_name == 'uv_pq_te_from_tm_with_far') then
            allocate(UVPQTEFromTMWithFarModel::model)
        else
            call assert(.false., 'unknown model name '//model_name)
        endif
    end function construct_calculation_model

end module calculation_models