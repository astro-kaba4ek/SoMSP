module far_context_module
    use regime
    use scattering_context_module
    use spheroidal_context_module
    use legendre_functions

    implicit none
    ! for transition between spheroidal basis and spherical nonorthogonal ('far') basis
    type :: FarContext
        integer :: m, spheroidal_lnum, spherical_lnum
        ! [1..spheroidal_lnum,1..spherical_lnum]
        complex(knd), allocatable, dimension(:, :) :: connecting_matrix
        integer :: state
    contains
        procedure :: initialize => initialize_far_context
        procedure :: reset => reset_far_context
        final :: delete_far_context
    end type FarContext

contains
    subroutine reset_far_context(this)
        class(FarContext), intent(out) :: this

        this%state = 0
    end subroutine reset_far_context

    subroutine initialize_far_context(this, m, spheroidal_lnum, spherical_lnum, outside_layer)
        class(FarContext), intent(inout) :: this
        integer, intent(in) :: m, spheroidal_lnum, spherical_lnum
        type(SpheroidalCalculation), intent(in) :: outside_layer

        real(knd) :: start, finish

        if (.not. allocated(this%connecting_matrix) .or. &
            any([m, spheroidal_lnum, spherical_lnum] /= [this%m, this%spheroidal_lnum, this%spherical_lnum])) then
            this%state = 0
        endif

        if (this%state == 1) then
            return
        endif

        call cpu_time(start)

        if (allocated(this%connecting_matrix) .and. any([spheroidal_lnum, spherical_lnum] /= shape(this%connecting_matrix))) then
            deallocate(this%connecting_matrix)
        endif

        if (.not. allocated(this%connecting_matrix)) then
            allocate(this%connecting_matrix(spheroidal_lnum, spherical_lnum))
        end if

        this%m = m
        this%spheroidal_lnum = spheroidal_lnum
        this%spherical_lnum = spherical_lnum

        call set_connecting_matrix(outside_layer, spheroidal_lnum, spherical_lnum, this%connecting_matrix)

        this%state = 1

        call cpu_time(finish)

        call log_time('context far', finish - start)

    end subroutine initialize_far_context

    subroutine set_connecting_matrix(layer, nsize, lsize, res)
        type(SpheroidalCalculation), intent(in) :: layer
        integer, intent(in) :: lsize, nsize
        complex(knd), intent(out) :: res(nsize, lsize)

        real(knd) :: legendre_norm(lsize)
        integer :: i, j

        if (nsize > layer%lnum) then
            write(*,*) '{ERROR}: legendre coeffcient size is too small: nsize = ', nsize, 'lnum = ', layer%lnum
            call exit(1)
        end if

        legendre_norm = calculate_legendre_norm(layer%m, lsize)

        do i = 1, nsize
            do j = 1, lsize
                res(i, j) = cmplx(0q0, 1q0, knd) ** (j - i) * layer%legendre(j - 1, i) * legendre_norm(j)
            end do
        end do

    end subroutine set_connecting_matrix

    subroutine delete_far_context(this)
        type(FarContext), intent(inout) :: this

        if (allocated(this%connecting_matrix)) then
            deallocate(this%connecting_matrix)
        endif
    
    end subroutine delete_far_context
end module far_context_module