!  Contains subroutines for matrix operations.
!  All functions utilize complex matrices with the kind parameter
!  set in the regime module from regime.f90 as knd.
module matrix
    use regime
    implicit none

    ! Matrix inversion in O(n^3) time
    ! needs more testing
    !private quick_inverse_matrix
    ! private inverse_matrix

contains

    !  Swaps rows i and j of the matrix A nxn
    !  Used in LU_gaussian
    subroutine swap_rows(A, n, i, j)
        complex(knd) A(n, n), v
        integer :: n, i, j, k

        do k = 1, n
            v = A(i, k)
            A(i, k) = A(j, k)
            A(j, k) = v
        enddo

    end subroutine swap_rows

    subroutine swap_columns(A, n, i, j)
        complex(knd) A(n, n), v
        integer :: n, i, j, k

        do k = 1, n
            v = A(k, i)
            A(k, i) = A(k, j)
            A(k, j) = v
        enddo

    end subroutine swap_columns

    !  Performs LU decomposition of the matrix A nxn to lower L and upper U
    !  with P as a row permutation matrix
    !  If the matrix is degenerate prints message to the console
    subroutine LU_gaussian(A, n, L, U, P)
        integer :: n, i, j, k
        complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), v, x(n, n), y(n, n)
        !  if abs(a(i,i)) < EPS it is considered to be =0 and
        !  hence the matrix is degenerate
        real(knd), parameter :: EPS = 1q-32
        complex(knd) :: v1, v2

        U = A
        L = 0
        call get_identity_matrix(P, n)

        do k = 1, n
            !            write(*,*) 'check-1: k = ',k
            !            write(*,*) 'P'
            !            do i = 1, n
            !                write(*,*) p(i,:)
            !            end do
            !            write(*,*) 'A'
            !            do i = 1, n
            !                write(*,*) A(i,:)
            !            end do
            !            x = matmul(P, A)
            !            write(*,*) 'matmul(P, A)'
            !            do i = 1, n
            !                write(*,*) x(i,:)
            !            end do
            !            write(*,*) 'L'
            !            do i = 1, n
            !                write(*,*) L(i,:)
            !            end do
            !            write(*,*) 'U'
            !            do i = 1, n
            !                write(*,*) U(i,:)
            !            end do
            !            write(*,*) 'matmul(L, U)'
            !            y = matmul(L, U)
            !            do i = 1, n
            !                write(*,*) y(i,:)
            !            end do
            L(k, k) = cmplx(1q0, 0q0, knd)

            if (abs(u(k, k)) < EPS) then
                do i = k + 1, n
                    if (abs(u(i, k)) > EPS) then
!                        write(*,*) 'swapping'
                        call swap_rows(U, n, i, k)
                        call swap_columns(P, n, i, k)
                        exit
                    endif
                enddo
                if (i > n) then
                    write(*, *) 'matrix is degenerate, u(k,k) = ', u(k, k)
                endif
            endif
!            write(*,*) 'start substracting'
            do i = k + 1, n
!                write(*,*) 'i = ', i
                L(i, k) = u(i, k) / u(k, k)
                v1 = u(i,k)
                v2 = u(k,k)
                !write(*, *) 'i = ', i, 'k = ', k, 'l(i,k)=', l(i, k)
                do j = 1, n
                    if (i == 53 .and. j == 1) then
!                    write(*,*) 'j = ', j, 'uij = ', u(i,j), 'lik = ', l(i,k), 'uik = ', u(i,k), 'ukk = ', u(k,k), 'lu = ', &
!                            u(k,k) * l(i,k), 'prod = ', L(i, k) * u(k, j),&
!                            'res = ', u(i, j) - L(i, k) * u(k, j)
!                        write(*,*) 'uik = ', u(i,k), 'ukj = ', u(k,j), 'pr = ', u(i,k)*u(k,j), 'ukk = ', u(k,k),&
!                                'res = ', u(i,k)*u(k,j) / u(k,k)

                            end if
!                    u(i, j) = u(i, j) - u(k, j) * (u(i, k) / u(k,k))
                    u(i, j) = (u(i, j) * v2 - u(k, j) * v1) / v2
!                    u(i, j) = u(i, j) - L(i, k) * u(k, j)
                enddo
!                u(i,k) = 0
            enddo
!            write(*,*)
        enddo

!        v1 = cmplx(1.65q36, 20q0, knd)
!        v2 = 1q0 / v1
!        write(*,*) v1, v2, v1 * v2

    end subroutine LU_gaussian

    !   Solves the linear system Ax = b, A:nxn, b:n
    !   result in the last parameter res
    subroutine solve_system(A, n, b, res)
        integer :: n, i, j, k
        complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), b(n), res(n), res_low(n)

        call LU_gaussian(A, n, L, U, P)
        b = matmul(P, b)

        call solve_lower_system(L, n, b, res_low)
        call solve_upper_system(U, n, res_low, res)

    end subroutine solve_system

    !  Solves a linear system with a lower triangular matrix
    subroutine solve_lower_system(A, n, b, res)
        integer :: n, i, j
        complex(knd) :: A(n, n), b(n), res(n)

        res = b

        do i = 1, n
            do j = 1, i - 1
                res(i) = res(i) - res(j) * A(i, j)
            enddo
            res(i) = res(i) / A(i, i)
        enddo

    end subroutine solve_lower_system

    !  Solves a linear system with an upper triangular matrix
    subroutine solve_upper_system(A, n, b, res)
        integer :: n, i, j
        complex(knd) :: A(n, n), b(n), res(n)

        res = b

        do i = n, 1, -1
            do j = i + 1, n
                res(i) = res(i) - res(j) * A(i, j)
            enddo
            res(i) = res(i) / A(i, i)
        enddo

    end subroutine solve_upper_system

    ! Calculates A^(-1) in O(n^4) operations where A:nxn
    subroutine inverse_matrix(A, n, res, sw)
        integer :: n, i, j, k
        logical, optional :: sw
        logical :: swap
        complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), res(n, n), res_low(n), &
                x(n, n), y(n, n)

        swap = .false.
        if (present(sw)) then
            swap = sw
        end if
        if (swap) then
            do i = 1, n, 2
                call swap_rows(A, n, i, i + 1)
            end do
        end if
!        A = A * cmplx(1q0, 1q0, knd)
        call LU_gaussian(A, n, L, U, P)
!                write(*,*) 'check0:'
!                write(*,*) 'P'
!                do i = 1, n
!                    write(*,*) p(i,:)
!                end do
!                write(*,*) 'A'
!                do i = 1, n
!                    write(*,*) A(i,:)
!                end do
        !        x = matmul(P, A)
        !        write(*,*) 'matmul(P, A)'
        !        do i = 1, n
        !            write(*,*) x(i,:)
        !        end do
!                write(*,*) 'L'
!                do i = 1, n
!                    write(*,*) L(i,:)
!                end do
!                write(*,*) 'U'
!                do i = 1, n
!                    write(*,*) U(i,:)
!                end do
!                write(*,*) 'p*matmul(L, U)'
!                y = matmul(P, matmul(L, U))
!                do i = 1, n
!                    write(*,*) y(i,:)
!                end do

        do k = 1, n

            call solve_lower_system(L, n, P(k, :), res_low)
!            write(*, *) 'check1:'
!            write(*, *) P(k, :)
!            write(*, *) matmul(L, res_low)
            call solve_upper_system(U, n, res_low, res(:, k))
!            write(*, *) 'check2:'
!            write(*, *) res_low
!            write(*, *) matmul(U, res(:, k))

        enddo
        if (swap) then
            do i = 1, n, 2
                call swap_rows(A, n, i, i + 1)
                call swap_columns(res, n, i, i + 1)
            end do
        end if
!        write(*, *) 'end check:'
!        x = matmul(A, res)
!        do i = 1, n
!            write(*, *) x(i, :)
!        end do
!        res = res * cmplx(1q0, 1q0, knd)

    end subroutine inverse_matrix

    !  Calculates the inverse matrix and determinant in O(n^3) operations where
    !  A:nxn. Does not perform correctly if A(i,i)=0 is encountered.
    subroutine quick_inverse_matrix(A, n, res)
        integer :: n, i, j, k
        complex(knd) :: A(n, n), L(n, n), U(n, n), P(n, n), res(n, n), res_low(n), det, mul
        !  if abs(a(i,i)) < EPS it is considered to be =0 and
        !  hence the matrix is degenerate
        real(knd), parameter :: EPS = 1q-32

        det = 1q0

        res = A

        do k = 1, n
!            write(*,*) 'k = ', k, 'res(k,k) = ', res(k,k)
            if (abs(res(k, k)) < EPS) then
                write(*, *) 'k = ', k, 'res(k,k) = ', res(k, k)
                !exit
                cycle
            end if

            det = det * res(k, k)
!            write(*,*) 'det = ', det
            mul = 1.0_knd / res(k, k)
            do i = 1, n
                if (i /= k) then
                    res(i, k) = -res(i, k) * mul
                    !                    res(k, i) = res(k, i) / res(k, k)
                end if
            end do

            do i = 1, n
                do j = 1, n
                    if (i /= k .and. j /= k) then
                        res(i, j) = res(i, j) + res(k, j) * res(i, k)
                    end if
                end do
            end do

            do i = 1, n
                if (i /= k) then
                    res(k, i) = res(k, i) * mul
                end if
            end do

            res(k, k) = mul
        enddo

    end subroutine quick_inverse_matrix

    ! Performs A := AxD, where A is a given matrix nxn, D - diagonal matrix
    ! which diagonal elements are given in diag_matr:1xn
    subroutine multiply_by_diag_right(A, n, diag_matr)
        complex(knd) A(n, n), diag_matr(n)
        integer :: n, i, j

        do i = 1, n
            do j = 1, n
                A(i, j) = A(i, j) * diag_matr(j)
            enddo
        enddo
    end subroutine multiply_by_diag_right

    ! Performs A := DxA, where A is a given matrix nxn, D - diagonal matrix
    ! which diagonal elements are given in giag_matr:1xn
    subroutine multiply_by_diag_left(A, n, diag_matr)
        complex(knd) A(n, n), diag_matr(n)
        integer :: n, i, j

        do i = 1, n
            do j = 1, n
                A(i, j) = diag_matr(i) * A(i, j)
            enddo
        enddo
    end subroutine multiply_by_diag_left

    !
    subroutine check_symmetric(A, n)
        integer :: n, i, j
        complex(knd) :: A(n, n)
        real(knd) :: rel_diff, abs_diff

        rel_diff = 0q0
        abs_diff = 0q0
        do i = 2, n
            !			do j = i + 1, n
            !				if (cqabs(A(i,j) + A(j,i)) > 1q-31) then
            !					if (A(i, j) /= A(j, i)) then
            !						write(*,*) 'matr', A(i,j), A(j,i)
            !					endif
            !					abs_diff = max(abs_diff, cqabs(A(i, j) - A(j, i)))
            !					rel_diff = max(rel_diff, cqabs(A(i, j) - A(j, i)) / cqabs(A(i,j) + A(j,i)))
            !				endif
            !			enddo
            abs_diff = max(abs_diff, abs(A(i, 1) - A(1, i)))
            rel_diff = max(rel_diff, abs(A(i, 1) - A(1, i)) / abs(A(i, 1) + A(1, i)))
        enddo
        write(*, *) 'maximum symmetric absolute difference in the first row = ', abs_diff
        write(*, *) 'maximum symmetric relative difference in the first row = ', rel_diff
    end subroutine check_symmetric

    !
    subroutine check_symmetric_result(A, n, abs_diff, rel_diff)
        integer :: n, i, j
        complex(knd), intent(in) :: A(n, n)
        real(knd), intent(out) :: rel_diff, abs_diff

        rel_diff = 0q0
        abs_diff = 0q0
        do i = 2, n
            abs_diff = max(abs_diff, abs(A(i, 1) - A(1, i)))
            rel_diff = max(rel_diff, abs(A(i, 1) - A(1, i)) / abs(A(i, 1) + A(1, i)))
        enddo
    end subroutine check_symmetric_result

    !  A := E_n
    subroutine get_identity_matrix(A, n)
        integer :: n, i
        complex(knd) :: A(n, n)
        A = 0
        do i = 1, n
            A(i, i) = cmplx(1q0, 0q0, knd)
        end do
    end subroutine get_identity_matrix

    ! A(i,i) := R(i), A(i,j) := 0, i != j
    subroutine get_full_matrix_from_diag(R, n, A)
        integer :: n, i
        complex(knd) :: R(n), A(n, n)

        A = 0
        do i = 1, n
            A(i, i) = R(i)
        end do
    end subroutine get_full_matrix_from_diag

    ! calculates the L_1 norm of a matrix
    real(knd) function norm_l1(A, n)
        integer :: n, i
        complex(knd) :: A(n, n)

        norm_l1 = 0q0
        do i = 1, n
            norm_l1 = max(norm_l1, sum(abs(A(:, i))))
        end do
    end function norm_l1

    ! calculates the L_\infty norm of a matrix
    real(knd) function norm_linf(A, n)
        integer :: n, i
        complex(knd) :: A(n, n)

        norm_linf = 0q0
        do i = 1, n
            norm_linf = max(norm_linf, sum(abs(A(i, :))))
        end do
    end function norm_linf

    real(knd) function max_eigen_triangle(A, n)
        integer :: n, i
        complex(knd) :: A(n, n)

        max_eigen_triangle = 0q0

        do i = 1, n
            max_eigen_triangle = max(max_eigen_triangle, abs(A(i, i)))
        end do
    end function max_eigen_triangle

    real(knd) function norm_l2(A, n)
        integer :: n, i
        complex(knd) :: A(n, n), B(n, n)
        complex(knd) :: L(n, n), U(n, n), P(n, n)

        B = transpose(conjg(A))
        B = matmul(B, A)

        call LU_gaussian(B, n, L, U, P)

        norm_l2 = sqrt(max_eigen_triangle(L, n) * max_eigen_triangle(U, n))
    end function norm_l2

    real(knd) function norm_l2_2(A, n)
        integer :: n, i, j
        complex(knd) :: A(n, n)

        norm_l2_2 = 0

        do i = 1, n
            do j = 1, n
                norm_l2_2 = norm_l2_2 + abs(A(i, j))**2
            end do
        end do

        norm_l2_2 = sqrt(norm_l2_2)

    end function norm_l2_2

    !real(knd) function max_eigen(A, n)
    !    EXTERNAL         ZGEEV
    !    integer :: n, i, j
    !    complex(knd) :: A(n,n)

    !    INTEGER :: LDA, LDVL, LDVR
    !    INTEGER, PARAMETER :: LWMAX = 1000
    !    INTEGER          INFO, LWORK
    !    DOUBLE PRECISION RWORK( 2*N )
    !    COMPLEX*16       A8( n, N ), VL( n, N ), VR( n, N ), W( N ), WORK( LWMAX )
    !    A8 = A
    !    LDA = N
    !    LDVL = N
    !    LDVR = N

    !    LWORK = -1
    !    CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
    !    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    !    CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )

    !    IF( INFO.GT.0 ) THEN
    !        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !        return
    !    END IF

    !    max_eigen = 0q0
    !    do i = 1, n
    !        max_eigen = max(max_eigen, abs(w(i)))
    !    end do

    !end function max_eigen

    !real(knd) function norm_l2_lapack(A, n)
    !    integer :: n
    !    complex(knd) :: A(n,n), B(n,n)

    !    B = transpose(conjg(A))
    !    B = matmul(B, A)

    !    norm_l2_lapack = qsqrt(max_eigen(B, n))
    !end function norm_l2_lapack

end module matrix
