!  Contains various functions for conversion of geometric properties of a spheroidal scatterer
module geometry
    use regime
    use constants
    implicit none
contains
    !  for the given rv and ab returns a - the major semiaxis
    real(knd) function get_geom_a(f, rv, ab) result(a)
        integer :: f
        real(knd) :: rv, ab

        if (f == 1) then
            a = rv * (ab ** (2q0 / 3q0))
        else
            a = rv * (ab ** (1q0 / 3q0))
        endif
    end function get_geom_a

    ! !  for the given rv and ab returns b - the minor semiaxis
    ! real(knd) function get_b(f, rv, ab) result(b)
    !     integer :: f
    !     real(knd) :: rv, ab

    !     if (f == 1) then
    !         b = rv / ab ** (1q0 / 3q0)
    !     else
    !         b = rv / ab ** (2q0 / 3q0)
    !     endif
    ! end function get_b

    !  Empiric function, evaluates the required size of matrices for the accuracy 1e-6
    integer function size_of_matrices(f, xv, ab, ri)
        integer, intent(in) :: f
        real(knd), intent(in) :: xv, ab
        complex(knd), intent(in) :: ri

        size_of_matrices = aint(get_geom_a(f, xv, ab) * max(1.2, real(ri))) + 8
        if (mod(size_of_matrices, 2) == 1) then
            size_of_matrices = size_of_matrices + 1
        end if

    end function size_of_matrices

    !  Expected maximum value of m for iteration
    !  Calculation stops if the next summand is less than a fraction
    !  MIN_M_RATIO from module constatnt of an already calculated
    !  value, but all considered m are less that this value
    integer function get_max_m(f, xv, ab)
        integer, intent(in) :: f
        real(knd), intent(in) :: xv, ab

        get_max_m = aint(get_geom_a(f, xv, ab)) * 2 + 10
    end function get_max_m

end module geometry