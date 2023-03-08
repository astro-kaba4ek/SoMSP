!  contains general constants and functions for obtaining parameter values
!  like size of matrices and number of legendre coefficients
module constants
    use regime
    use mpi_f08

    implicit none

    ! MPI constants
    integer, public :: Err, Rank, Size_mpi
    TYPE(MPI_Status) :: Status

    real(knd), parameter :: PI = 4.0_knd * atan(1.0_knd)
    real(knd), parameter :: BASE_ACCURACY = 1e-8_knd
    real(knd), parameter :: MIN_M_RATIO = 1e-8_knd
    real(knd), parameter :: DEFAULT_PAIR_INTEGRAL_ACCURACY = 1e-24_knd;
    real(4), parameter :: REVERSE_MATRIX_SIZE_RATIO = 1.0 !116.0 / 102.0
    integer, parameter :: MAXIMUM_MATRIX_SIZE = 400
    integer, parameter :: MAXIMUM_D_COEFF_NUMBER = 100000
    integer, parameter :: BUCKET_SIZE = 20
    integer, parameter :: GLOBAL_BUCKET = 10000000
    real(knd), parameter :: THETA_ADD = 1.0e-8_knd
    complex(knd), parameter :: IDEG(0:3) = (/ &
            cmplx(1.0_knd, 0.0_knd, knd), cmplx(0.0_knd, 1.0_knd, knd), &
            cmplx(-1.0_knd, 0.0_knd, knd), cmplx(0.0_knd, -1.0_knd, knd)/)
    complex(knd), parameter :: NEGIDEG(0:3) = (/ &
            cmplx(1.0_knd, 0.0_knd, knd), cmplx(0.0_knd, -1.0_knd, knd), &
                    cmplx(-1.0_knd, 0.0_knd, knd), cmplx(0.0_knd, 1.0_knd, knd)/)

    enum, bind(c)
        ! Spheroidal type
        enumerator :: OBLATE = -1, PROLATE = 1
        !  levels of logs
        !  DO NOT CHANGE THE ORDER! PARAMETER ARRAYS FOR LOGS ARE DEFINED AS
        !  ARRAY(ERROR:DETAIL)
        enumerator :: ERROR, WARNING, DEBUG, INFO, DETAIL
        ! Scattering modes
        enumerator TE, TM, TETM
        ! Bases
        enumerator SPHEROIDAL_BASIS, FAR_BASIS, BARBER_BASIS, MISHCH_BASIS
        ! Potentials
        enumerator UV, PQ, UVPQ
        ! Types of scattering factors
        enumerator :: qfactors, cfactors, normalized_cfactors
        enumerator :: no_absorbtion, consequential
    end enum

    integer, parameter :: SCAT_MATR_FD = 14
    integer, parameter :: LOG_FD = 16
    character(128), parameter :: LOG_FILENAME = 'scattering.log'

#ifdef NEED_LOG
    logical, parameter :: LOG_INFO = .true.
#else
    logical, parameter :: LOG_INFO = .false.
#endif

#ifdef NEED_LOG_TIME
    logical, parameter :: LOG_TIMES = .true.
#else
    logical, parameter :: LOG_TIMES = .false.
#endif

#ifdef NEED_LOG_TMATRIX
    logical, parameter :: LOG_MATRICES = .true.
#else
    logical, parameter :: LOG_MATRICES = .false.
#endif

#ifdef NEED_LOG_AMPLITUDE
    logical, parameter :: LOG_AMPLITUDE_MATRIX = .true.
#else
    logical, parameter :: LOG_AMPLITUDE_MATRIX = .false.
#endif

    logical, parameter :: LOG_BLOCKS = .false.
contains

    integer function get_full_matrix_size(matrix_size)
        integer, intent(in) :: matrix_size
        get_full_matrix_size = aint(REVERSE_MATRIX_SIZE_RATIO * matrix_size)
    end function get_full_matrix_size

    integer function get_full_function_size(matrix_size, spherical_lnum)
        integer, intent(in) :: matrix_size, spherical_lnum
        get_full_function_size = max(aint(REVERSE_MATRIX_SIZE_RATIO * matrix_size), 1.0 * spherical_lnum)
    end function get_full_function_size

end module constants