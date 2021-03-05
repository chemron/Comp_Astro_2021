module init
    implicit none
    real , parameter :: pi = 4.* atan (1.0)

    contains
    subroutine setup(x, v, m, h, rho, u, P, c, n_max, n)

        integer, intent(in) :: n_max
        integer, intent(out) :: n
        real, parameter :: rho_0 = 1.0
        real, parameter :: x_min = 0.0, x_max = 1.0
        real :: dx
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
        n = 100

        ! setup position
        dx = (x_max - x_min)/(n - 1)
        do i=1, n
            x(i) = x_min + dx*(i - 1)
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
            m(i) = rho_0*dx
        enddo

        ! setup smoothing length
        do i = 1, n
            h(i) = 1.2 * dx
        enddo

        ! setup density

    end subroutine setup

    subroutine set_ghosts(x, v, m, h, rho, u, P, c, n_max, n_ghosts, n)

        integer, intent(in) :: n_max, n_ghosts, n
        real :: dx
        integer :: i
        real, parameter :: rho_0 = 1.0
        real, intent(inout) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)

        ! setup position
        dx = (x(100) - x(1))/(n - 1)
        do i=101, 101 + n_ghosts/2 - 1
            x(i) = x(i-1) + dx
            c(i) = 1
            v(i) = 0
            m(i) = rho_0*dx
            h(i) = 1.2 * dx
        enddo

        do i= 1, n_ghosts/2
            x(i + 101 + n_ghosts/2 - 1) = -dx*i
            c(i + 101 + n_ghosts/2 - 1) = 1
            v(i + 101 + n_ghosts/2 - 1) = 0
            m(i + 101 + n_ghosts/2 - 1) = rho_0*dx
            h(i + 101 + n_ghosts/2 - 1) = 1.2 * dx
        enddo

    end subroutine set_ghosts


    subroutine output(x, v, m, h, rho, u, P, c, n_max, n, n_ghosts)
        integer, intent(in) :: n_max, n, n_ghosts
        real, intent(in) ::  x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max)
        integer :: lu = 1, i
        real :: t = 0

        open(lu , file='output.out', status='replace', action='write')
        write(lu,*) '# x, v, m, h, rho, u, P, c'
        write(lu,*) t
        do i=1,n + n_ghosts
            write(lu,*) x(i), v(i), m(i), h(i), rho(i), u(i), P(i), c(i)
        enddo

    end subroutine output

end module init
