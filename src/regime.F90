!  contains kind parameters for general calculations
module regime
#ifdef USE_DOUBLE_PRECISION
    integer, parameter :: knd = 8
#else
    integer, parameter :: knd = 16
#endif
end module regime