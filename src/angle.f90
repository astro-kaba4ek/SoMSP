module angle
    use regime
    use constants
    implicit none
    private
    !  structure intended to optimize cos and sin calculation
    type, public :: AngleType
        real(knd) :: value
        real(knd) :: angle_cos
        real(knd) :: angle_sin
    contains
        private
        procedure, public :: set => set_angle
        procedure, public :: in_degrees
    end type AngleType

    interface AngleType
        procedure :: construct_angle_type
    end interface AngleType

contains

    subroutine set_angle(this, value)
        class(AngleType), intent(out) :: this
        real(knd), intent(in) :: value

        this%value = value
        this%angle_cos = cos(value)
        this%angle_sin = sin(value)
    end subroutine set_angle

    function construct_angle_type(value) result(res)
        real(knd), intent(in) :: value
        type(AngleType) :: res

        call res%set(value)
        
    end function construct_angle_type

    real(knd) function in_degrees(this)
        class(AngleType), intent(in) :: this

        in_degrees = this%value * 180.0_knd / PI
    end function in_degrees

end module angle