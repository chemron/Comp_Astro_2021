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
        integer :: i
        real, parameter :: rho_0 = 1.0
        real, intent(inout) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
        real :: width

        ! get width
        width = abs(x(n) - x(1))

        ! setup position
        do i=1, n_ghosts/2
            x(n + i) = x(i + 1) + width 
            c(n + i) = c(i + 1)
            v(n + i) = v(i + 1)
            m(n + i) = m(i + 1)
            h(n + i) = h(i + 1)
        enddo

        do i= 1, n_ghosts/2
            x(i + n + n_ghosts/2) = x(n - i) - width
            c(i + n + n_ghosts/2) = c(n - i)
            v(i + n + n_ghosts/2) = v(n - i)
            m(i + n + n_ghosts/2) = m(n - i)
            h(i + n + n_ghosts/2) = h(n - i)
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
