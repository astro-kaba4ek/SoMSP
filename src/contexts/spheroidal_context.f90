module spheroidal_context_module
    use regime
    use constants
    use integrals
    use matrix
    use scattering_context_module
    use logging
    use legendre_functions

    implicit none

    type :: SpheroidalContext
        ! general
        ! fnum - number of function values to calculate
        ! lnum - number of field/potential expansion 
        ! must be fnum >= lnum
        integer :: m, fnum, lnum, nol
        ! functions
        ! layers(0,j) = outside_layer for layer j (with ref_ind related to the outside and c, ksi to this layer's size
        ! layers(1,j) - inside_layer
        ! [0..1][1..nol]
        type(SpheroidalCalculation), allocatable :: layers(:,:)
        integer, private :: maxd
        ! integrals [1..lnum, 1..lnum, 1..nol]
        ! PQ
        complex(knd), allocatable, dimension(:, :, :) :: Delta
        ! UV
        complex(knd), allocatable, dimension(:,:,:) ::  Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, Pi1, Pi3
        ! state: 0 - not initialized, 1 - only functions are calculated, 2 - for PQ, 3 - for UV
        integer :: state
    contains
        procedure, private :: check_spheroidal_functions, check_base_matrices
        procedure :: initialize => initialize_spheroidal_context
        procedure :: reset => reset_spheroidal_context
        final :: delete_spheroidal_context
    end type SpheroidalContext

contains
    subroutine reset_spheroidal_context(this)
        class(SpheroidalContext), intent(out) :: this

        this%state = 0
    end subroutine reset_spheroidal_context

    subroutine initialize_spheroidal_context(this, m, fnum, lnum, scattering_context, target_state)
        class(SpheroidalContext), intent(inout) :: this
        ! expansion parameters
        integer, intent(in) :: m, fnum, lnum
        ! physical scattering state for functions arguments
        type(ScatteringContext), intent(in) :: scattering_context
        ! what needs to be calculated: 
        !   1 - only functions, 
        !   2 - functions and delta integral matrix (for pq solution)
        !   3 - functions and all integral matrices (for uv solution)
        integer, intent(in) :: target_state

        if (target_state < 1 .or. target_state > 3) then
            write(LOG_FD, *) '{ERROR} invalid target_state = ', target_state
            call exit(1)
        endif

        if (fnum < lnum) then
            write(LOG_FD, *) '{ERROR} fnum = ', fnum, ' is less than lnum = ', lnum
            call exit(1)
        endif

        ! allocated check is used to see whether the state variable is valid, if it is not, it is equivalent to it being 0
        if (.not. allocated(this%layers) .or. m /= this%m .or. &
            any([fnum, lnum, scattering_context%scatterer%number_of_layers] > [this%fnum, this%lnum, this%nol])) then
            this%state = 0
        endif

        if (this%state >= target_state) return

        ! if (LOG_INFO) write(LOG_FD,*) 'initialize spheroidal context with target state =', target_state

        if (this%state == 0) then
            this%m = m
            this%fnum = fnum
            this%lnum = lnum
            this%nol = scattering_context%scatterer%number_of_layers

            call check_spheroidal_functions(this, scattering_context)

            this%state = 1
        endif

        if (target_state > 1) then
            call check_base_matrices(this, target_state)
        endif

    end subroutine initialize_spheroidal_context

    subroutine check_spheroidal_functions(this, scattering_context)
        class(SpheroidalContext), intent(inout) :: this
        type(ScatteringContext), intent(in) :: scattering_context

        integer :: j
        real(knd) :: start, finish

        if (this%state /= 0) return

        if (allocated(this%layers) .and. (size(this%layers(0,:)) /= this%nol)) then
            deallocate(this%layers)
        end if
        if (.not. allocated(this%layers)) then
            allocate(this%layers(0:1, 1:this%nol))
        end if

        call cpu_time(start)

        ! layer(0, j) outside layer of the j layer of the particle
        ! layer(1, j) - inside layer of j layer
        ! j = 1 - mantle
        ! j = 2 - core
        do j = 1, this%nol
            call this%layers(0,j)%calculate( &
                this%m, &
                this%fnum, &
                scattering_context%scatterer%c0(j) * scattering_context%calculation_point%get_refractive_index(j - 1), &
                scattering_context%scatterer%ksi(j), &
                1, &
                (/ scattering_context%directions%alpha%angle_cos /), &
                scattering_context%scatterer%spheroidal_type &
            )
            call this%layers(1,j)%calculate( &
                this%m, &
                this%fnum, &
                scattering_context%scatterer%c0(j) * scattering_context%calculation_point%get_refractive_index(j), &
                scattering_context%scatterer%ksi(j), &
                1, &
                (/ scattering_context%directions%alpha%angle_cos /), &
                scattering_context%scatterer%spheroidal_type &
            )
        end do

        this%state = 1

        call cpu_time(finish)

        call log_time('context spheroidal functions', finish - start)

    end subroutine check_spheroidal_functions

    subroutine calculate_base_matrices(layer0, layer1, matrix_size, &
        Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, mult_coef)
        type(SpheroidalCalculation), intent(in) :: layer0, layer1
        integer, intent(in) :: matrix_size
        complex(knd), intent(out) :: Delta(matrix_size, matrix_size), &
                Q01(matrix_size, matrix_size), Q11(matrix_size, matrix_size), Q01Q11(matrix_size, matrix_size), &
                Kappa(matrix_size, matrix_size), Gamma11(matrix_size, matrix_size), Epsilon(matrix_size, matrix_size)

        complex(knd), allocatable, dimension(:, :) :: tmp, result, identity
        real(knd) :: ksi
        integer :: m, full_size, i, j, k
        real(knd), intent(in) :: mult_coef(:)
        real(knd) :: start, finish

        if (layer0%m /= layer1%m) then
            write(*,*) 'different m in layers!'
            return
        endif

        m = layer1%m
        ksi = layer0%ksi

    !        call calculate_kappa(layer1, layer1, Kappa, matrix_size)
        call double_dep_integral(Kappa, m, layer1%legendre, layer1%legendre, mult_coef, kappa_c_lower, kappa_c_upper)
    !        call calculate_gamma(layer1, layer1, Gamma11, matrix_size)
        call double_dep_integral(Gamma11, m, layer1%legendre, layer1%legendre, mult_coef, gamma_c_lower, gamma_c_upper)
    !        call calculate_epsilon(layer1, layer1, Epsilon, matrix_size)
        call triple_dep_integral(Epsilon, m, layer1%legendre, layer1%legendre, mult_coef, epsilon_c_lower, epsilon_c_middle, &
                epsilon_c_upper)

        full_size = min(get_full_matrix_size(matrix_size), min(layer0%lnum, layer1%lnum))

        allocate(tmp(full_size, full_size), identity(full_size, full_size), result(full_size, full_size))

        call get_identity_matrix(identity, full_size)
    !        call calculate_omega(layer1, layer1, tmp, full_size)
        call triple_dep_integral(tmp, m, layer1%legendre, layer1%legendre, mult_coef, omega_c_lower, omega_c_middle, &
                omega_c_upper)
                
        tmp = (ksi**2 - layer0%spheroidal_type) * identity + layer0%spheroidal_type * tmp
        call cpu_time(start)
        call quick_inverse_matrix(tmp, full_size, result)
        call cpu_time(finish)
        ! write(*,*) 'inverse time = ', finish - start
        Q11 = result(1:matrix_size, 1:matrix_size)

        identity = 0
    !        call calculate_delta(layer0, layer1, identity, full_size)
        call single_dep_integral(identity, m, layer0%legendre, layer1%legendre, mult_coef, delta_coef)
    !        Delta = 0
        Delta = identity(1:matrix_size, 1:matrix_size)
        call cpu_time(start)
        tmp = matmul(identity, result)
        ! do j = 1, full_size
        !     do k = 1, full_size
        !         do i = 1, matrix_size
        !             tmp(i,j) = tmp(i,j) + identity(i,k) * result(k,j)
        !         end do
        !     end do
        ! end do
        call cpu_time(finish)
        ! write(*,*) 'first matmul time = ', finish - start
        Q01 = tmp(1:matrix_size, 1:matrix_size)
        call cpu_time(start)
        result = matmul(tmp, result)
        ! do j = 1, matrix_size
        !     do k = 1, full_size
        !         do i = 1, matrix_size
        !             Q01Q11(i,j) = Q01Q11(i,j) + tmp(i,k) * result(k,j)
        !         end do
        !     end do
        ! end do
        call cpu_time(finish)
        ! write(*,*) 'second matmul time = ', finish - start
        Q01Q11 = result(1:matrix_size, 1:matrix_size)

        deallocate(tmp, identity, result)

    end subroutine calculate_base_matrices

    subroutine mult_matr_by_i(a, lnum)
        complex(knd), intent(inout) :: a(lnum, lnum)
        integer, intent(in) :: lnum

        integer :: i, j

        do i = 1, lnum
            do j = 2 - mod(i,2), lnum, 2
                a(i,j) = a(i,j) * IDEG(mod(j-i + 4 * lnum, 4))
            enddo
        enddo
    end subroutine mult_matr_by_i

    subroutine check_base_matrices(this, target_state)
        class(SpheroidalContext), intent(inout) :: this
        integer, intent(in) :: target_state

        real(knd), allocatable :: mult_coef(:)
        integer :: j, nol, lnum, md, i, k
        real(knd) :: start, finish
        real(knd), allocatable :: p_coef(:)

        if (this%state >= target_state) return

        nol = this%nol
        lnum = this%lnum

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate integrals'
        if (LOG_INFO) write(LOG_FD, *) '{INFO} calculate base matrices of sizes ', lnum, 'x', lnum, 'x', nol

        if (allocated(this%Delta) .and. any(shape(this%Delta) /= (/lnum, lnum, nol/))) then
            deallocate(this%Delta, this%Pi1, this%Pi3)
            if (allocated(this%Q01)) then
                deallocate(this%Q01, this%Q11, this%Q01Q11, this%Kappa, this%Gamma11, this%Epsilon)
            endif
        end if

        if (.not. allocated(this%Delta)) then
            allocate(this%Delta(lnum, lnum, nol), this%Pi1(lnum, lnum, nol), this%Pi3(lnum, lnum, nol))
        end if

        if (target_state == 3) then
            if (.not. allocated(this%Q01)) then
            allocate(this%Q01(lnum, lnum, nol), this%Q11(lnum, lnum, nol), &
                    this%Q01Q11(lnum, lnum, nol), this%Kappa(lnum, lnum, nol), &
                    this%Gamma11(lnum, lnum, nol), this%Epsilon(lnum, lnum, nol))
            endif
        endif

        call cpu_time(start)

        do j = 1, this%nol
            md = max(this%layers(0, j)%maxd, this%layers(1, j)%maxd)
            p_coef = 1.0_knd / calculate_legendre_coef(this%m, md + 1)
            call fill_common_multiplier(this%m, md, mult_coef)

            if (target_state == 3) then
                call calculate_base_matrices(this%layers(0, j), this%layers(1, j), lnum, &
                        this%Delta(:,:,j), this%Q01(:,:,j), this%Q11(:,:,j), this%Q01Q11(:,:,j), &
                        this%Kappa(:,:,j), this%Gamma11(:,:,j), this%Epsilon(:,:,j), mult_coef)
            endif

            if (this%state < 2) then
                if (j < nol) then
                    call single_dep_integral(this%Pi1(:,:,j), this%m, this%layers(1,j)%legendre, this%layers(0,j + 1)%legendre, &
                            p_coef, delta_coef)
                    call mult_matr_by_i(this%Pi1(:,:,j), lnum)
                    this%Pi3(:,:,j) = this%Pi1(:,:,j)
                    call multiply_by_diag_left(this%Pi1(:,:,j), lnum, this%layers(1,j)%r1(1:lnum))
                    call multiply_by_diag_right(this%Pi1(:,:,j), lnum, 1.0_knd / this%layers(0,j + 1)%r1(1:lnum))
                    call multiply_by_diag_left(this%Pi3(:,:,j), lnum, this%layers(1,j)%r3(1:lnum))
                    call multiply_by_diag_right(this%Pi3(:,:,j), lnum, 1.0_knd / this%layers(0,j + 1)%r3(1:lnum))
                else
                    call get_identity_matrix(this%Pi1(:,:,j), lnum)
                    this%Pi3(:,:,j) = this%Pi1(:,:,j)
                end if

                if (target_state == 2) then
                    call single_dep_integral(this%Delta(:,:,j), this%m, this%layers(0, j)%legendre, this%layers(1, j)%legendre, &
                        mult_coef, delta_coef)
                endif
            endif
        end do

        deallocate(mult_coef, p_coef)

        this%state = max(this%state, target_state) 

        call cpu_time(finish)

        call log_time('context base matrices', finish - start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate integrals'

    end subroutine check_base_matrices


    subroutine delete_spheroidal_context(this)
        type(SpheroidalContext), intent(inout) :: this

        if (allocated(this%layers)) then
            deallocate(this%layers)
        endif

        if (allocated(this%Delta)) then
            deallocate(this%Delta)
        endif

        if (allocated(this%Q01)) then
            deallocate(this%Q01, this%Q11, this%Q01Q11, this%Kappa, this%Gamma11, this%Epsilon, this%Pi1, this%Pi3)
        endif

    end subroutine delete_spheroidal_context

end module spheroidal_context_module