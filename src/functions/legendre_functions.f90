! Created by drakosha on 18.02.2021.

module legendre_functions

    use regime
    use constants

    use complex_prolate_swf
!    use complex_oblate_swf

    implicit none

    type, public :: LegendreCalculation
        integer :: m
        !  lnum = n - m + 1 >= 0
        integer :: lnum
        !  number of arguments
        integer :: narg
        !  array of arguments [1..narg]
        real(knd), allocatable, dimension(:) :: args
        ! functions and their derivatives for n - m = 0,.., lnum - 1, arrays [1..narg][1:lnum]
        real(knd), allocatable, dimension(:,:) :: pr, pdr

        ! c(n) = (2n+1) / 2 * (n-m)!/(n+m)!
        ! coef(i) = 1/c(i+m-1) = 1/norm^2
        real(knd), allocatable, dimension(:) :: coef
        real(knd), allocatable, dimension(:) :: norm

        logical :: calculated
    contains
        !  initializes, allocates arrays, does not calculate
        procedure :: set
        !  call to calculate all the values, can be long
        procedure :: calculate => calculate_legendre_functions !, calculate_q
        !  destructor, deallocate arrays
        final :: delete_calculation
    end type LegendreCalculation

contains
    function calculate_legendre_coef(m, lnum) result(coef)
        integer, intent(in) :: m, lnum

        real(knd) :: coef(lnum)

        integer :: i

        coef(1) = 1q0
        do i = 1, 2 * m
            coef(1) = coef(1) / i
        end do

        do i = 2, lnum
            coef(i) = coef(i - 1) * (i - 1) / (i - 1 + m * 2)
        end do

        do i = 1, lnum
            coef(i) = coef(i) * (2.0_knd * (i - 1 + m) + 1) / 2.0_knd
        end do

    end function calculate_legendre_coef

    function calculate_legendre_norm(m, lnum) result(norm)
        integer, intent(in) :: m, lnum

        real(knd) :: norm(lnum)

        norm = calculate_legendre_coef(m, lnum)
        norm = 1.0_knd / sqrt(norm)

    end function calculate_legendre_norm

    subroutine set(this, m, n, eta)
        class(LegendreCalculation) :: this
        integer, intent(in) :: m, n
        real(knd), intent(in) :: eta

        this%m = m
        this%lnum = max(n - m + 1, 10)
        this%narg = 1
        if (.not.allocated(this%args)) then
            allocate(this%args(1))
        end if
        this%args(1) = eta

        this%calculated = .false.
    end subroutine set

    subroutine calculate_legendre_functions(this)
        class(LegendreCalculation), intent(inout) :: this

        integer :: i, j, maxp
        real(vb_knd), allocatable, dimension(:) :: alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe, pnorm, pdnorm
        real(vb_knd), allocatable, dimension(:,:) :: pr, pdr
        real(vb_knd) :: vb_arg(this%narg)
        integer, allocatable, dimension(:) :: ipdnorm, ipnorm

        maxp = this%lnum
        if (allocated(this%pr)) then
            deallocate(this%pr, this%pdr, this%coef)
        endif
        allocate(this%pr(1:this%narg, 0:maxp), this%pdr(1:this%narg, 0:maxp))
        allocate(alpha(maxp),beta(maxp),gamma(maxp),coefa(maxp),coefb(maxp),coefc(maxp),coefd(maxp),coefe(maxp))
        allocate(pr(this%narg, maxp), pdr(this%narg, maxp), pnorm(this%narg), pdnorm(this%narg))
        allocate(ipdnorm(maxp),ipnorm(maxp))
        if (allocated(this%coef)) then
            deallocate(this%coef)
        end if
        if (allocated(this%norm)) then
            deallocate(this%norm)
        end if
        allocate(this%coef(maxp), this%norm(maxp))

        this%coef(1) = 1q0
        do i = 1, 2 * this%m
            this%coef(1) = this%coef(1) / i
        end do
        do i = 2, maxp
            this%coef(i) = this%coef(i - 1) * (i - 1) / (i - 1 + this%m * 2)
        end do

        do i = 1, maxp
            this%coef(i) = (this%coef(i) * (2q0 * (i - 1 + this%m) + 1) / 2q0)
        end do

        this%coef = (this%coef)
        this%norm = 1q0 / (this%coef)

        vb_arg = this%args
        call pleg (this%m, this%lnum - 2, maxp, precision(vb_arg(1)), range(vb_arg(1)) - 1, 0, 1, &
                vb_arg, this%narg,this%narg,pr,&
                pdr,pdnorm,ipdnorm,pnorm,ipnorm, &
                alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
        this%pr(:,1:) = pr
        this%pdr(:,1:) = pdr
!        write(*,*) 'm = ', this%m
!        write(*,*) 'pnorm = ', pnorm
!        write(*,*) 'pdnorm = ', pdnorm
!        write(*,*) 'pr1 = ', pr(1,1:)
!        write(*,*) 'pr1*pnorm = ', pr(1,1:)*pnorm(1)
!        write(*,*) 'pdr1 = ', pdr(1,1:)
!        write(*,*) 'pdr1*pdnorm = ', pdr(1,1:)*pdnorm(1)
!        write(*,*)
        this%pr(:,1:) = pr
        this%pdr(:,1:) = pdr
!        this%pdr(:,0) = pdnorm
!        this%pr(:,0) = pnorm
!        this%pr(:,0) = this%pr(:,0) * (10q0**ipnorm)
!        this%pdr(:,0) = this%pdr(:,0) * (10q0**ipdnorm)
        do i = 1, this%narg
!            this%pr(i,1) = this%pr(i,1) * this%pr(i,0)
!            this%pdr(i,1) = this%pdr(i,1) * this%pdr(i,0)
            do j = 3, this%lnum
                this%pr(i, j) = this%pr(i, j) * this%pr(i, j - 2)
                this%pdr(i, j) = this%pdr(i, j) * this%pdr(i, j - 2)
            end do
        end do

        this%pr = this%pr * pnorm(1) * (10q0**ipnorm(1))
        this%pdr = this%pdr * pdnorm(1) * (10q0**ipdnorm(1))
        if (this%m == 0) then
            this%pdr(:, 0) = 0
            this%pdr(:, 1) = 0
        end if
        !write(*,*) 'pr = ', this%pr(1,1:)
        !write(*,*) 'pdr = ', this%pdr(1,1:)
        !write(*,*) 'coef = ', this%coef

        deallocate(alpha)
        deallocate(beta)
        deallocate(gamma)
        deallocate(coefa)
        deallocate(coefb)
        deallocate(coefc)
        deallocate(coefd)
        deallocate(coefe)
        deallocate(pr, pdr, pnorm, pdnorm)
        deallocate(ipdnorm,ipnorm)
        this%calculated = .true.

    end subroutine calculate_legendre_functions

!    subroutine calculate_q(this)
!        class(LegendreCalculation), intent(inout) :: this
!
!        integer :: i, j, maxp
!        integer :: iqdl(2 * this%lnum),iql(2 * this%lnum), iqdml, iqml, itermpq
!        real*16 qml,qdml,termpq
!        real*16, allocatable, dimension(:) :: qdl,qdr,ql,qr
!
!        maxp = this%lnum
!        if (allocated(this%pr)) then
!            deallocate(this%pr, this%pdr, this%coef)
!        endif
!        allocate(this%pr(1:this%narg, 0:maxp), this%pdr(1:this%narg, 0:maxp))
!        allocate(this%coef(maxp))
!        allocate(qdl(2 * this%lnum),qdr(2 * this%lnum),ql(2 * this%lnum),qr(2 * this%lnum))
!
!        this%coef(1) = 1q0
!        do i = 1, 2 * this%m
!            this%coef(1) = this%coef(1) / i
!        end do
!        do i = 2, maxp
!            this%coef(i) = this%coef(i - 1) * (i - 1) / (i - 1 + this%m * 2)
!        end do
!
!        do i = 1, maxp
!            this%coef(i) = this%coef(i) * (2q0 * (i - 1 + this%m) + 1) / 2q0
!        end do
!
!        call qleg (this%m,this%lnum,this%lnum-2,this%lnum,real(this%args(1), 16),31,&
!                qdr,qdml,iqdml,qdl,&
!                iqdl,qr,qml,iqml,ql,iql,termpq,itermpq)
!
!        this%pr(1,1:this%lnum) = ql(1:this%lnum)
!        this%pdr(1,1:this%lnum) = qdl(1:this%lnum)
!        !this%pr(1,0) = qml * (10q0**iqml)
!        !this%pdr(1,0) = qdml * (10q0**iqdml)
!        !    do j = 1, this%lnum
!        !        this%pr(1, j) = qr(j) * this%pr(1, j - 1)
!        !        this%pdr(1, j) = qdr(j) * this%pdr(1, j - 1)
!        !    end do
!
!        !write(*,*) 'pr = ', this%pr(1,1:)
!        !write(*,*) 'pdr = ', this%pdr(1,1:)
!        !write(*,*) 'coef = ', this%coef
!
!        deallocate(qdl,qdr,ql,qr)
!
!        this%calculated = .true.
!
!    end subroutine calculate_q


    subroutine delete_calculation(this)
        type(LegendreCalculation), intent(inout) :: this

        if (allocated(this%pr)) then
            deallocate(this%pr, this%pdr, this%coef, this%norm)
        endif

        if (allocated(this%args)) then
            deallocate(this%args)
        end if

    end subroutine delete_calculation

end module legendre_functions