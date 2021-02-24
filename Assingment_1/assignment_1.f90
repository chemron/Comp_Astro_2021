module physics
    implicit none
    public :: setup

contains
    subroutine setup(x, v, m, h, rho, u, P, c, n_max, n)

        integer, intent(in) :: n_max
        integer, intent(out) :: n
        integer, parameter :: rho_0 = 1.0
        real, parameter :: x_min = 0.0, x_max = 1.0
        real , parameter :: pi = 4.* atan (1.0)
        real :: dx
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
        n = 100

        ! setup position
        dx = (x_max - x_min)/(n - 1)
        do i=0, n-1
            x(i + 1) = x_min + dx * i
        enddo

        ! setup sound speed
        ! TODO: make this proper
        do i = 1, n
            c(i) = 1
        enddo

        ! setup velocity
        do i = 1, n
            v(i) = 10.0 ** (-4) * c(i) * sin(2 * pi * (x(i)-x_min) / (x_max-x_min))
        enddo

        ! setup mass
        do i = 1, n
            m(i) = rho_0
        enddo

        do i = 1, n
            h(i) = 1.2 * dx
        enddo

    end subroutine setup


    subroutine output(x, v, m, h, rho, u, P, c, n_max, n)
        integer, intent(in) :: n_max, n
        real, intent(in) ::  x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max)
        integer :: lu = 1, i
        real :: t = 0

        open(lu , file='output.out', status='replace', action='write')
        write(lu,*) '# x, v, m, h, rho, u, P, c'
        write(lu,*) t
        do i=1,n
            write(lu,*) x(i), v(i), m(i), h(i), rho(i), u(i), P(i), c(i)
        enddo

    end subroutine output

end module physics


program assingment_1
    use physics

    implicit none
    integer, parameter :: n_max = 100
    integer :: n = 100
    real :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)
    print*, 'Hello World'


    call setup(x, v, m, h, rho, u, P, c, n_max, n)

    call output(x, v, m, h, rho, u, P, c, n_max, n)

end program assingment_1