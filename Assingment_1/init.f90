module init
    implicit none
    real , parameter :: pi = 4.* atan (1.0)

    contains
    subroutine setup(x, v, m, h, c_0, x_min, x_max, n_max, n)

        integer, intent(in) :: n_max
        integer, intent(out) :: n
        real, parameter :: rho_0 = 1.0
        real, intent(in) :: x_min, x_max
        real :: dx
        real, intent(in) :: c_0
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max)
        n = 100

        ! setup position
        ! put particles in the middle of cells with width dx
        dx = (x_max - x_min)/n
        do i=1, n
            x(i) = x_min + dx*(i - 0.5)
        enddo

        ! setup velocity
        do i = 1, n
            v(i) = 10.0 ** (-4) * c_0 * sin(2 * pi * (x(i)-x_min) / (x_max-x_min))
        enddo

        ! setup mass
        do i = 1, n
            m(i) = rho_0*dx
        enddo

        ! setup smoothing length
        do i = 1, n
            h(i) = 1.2 * dx
        enddo


    end subroutine setup

    subroutine isothermal_setup(x, v, m, h, x_min, x_max, n_max, n)

        integer, intent(in) :: n_max
        integer :: n_left, n_right
        integer, intent(out) :: n
        real, parameter :: rho_left = 1.0, rho_right = 0.1
        real, intent(inout) :: x_min, x_max
        real :: dx_left = 0.001, dx_right = 0.01
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max)

        ! calculate number of particles left and right of origin
        ! (round to nearest whole number)
        n_left = nint(abs(0 - x_min)/dx_left)
        n_right = nint(abs(x_max - 0)/dx_right)
        n = n_left + n_right

        ! left of origin
        do i=1, n_left
            ! position
            x(i) = x_min + dx_left*(i - 0.5)
            ! mass
            m(i) = rho_left*dx_left
            ! smoothing length
            h(i) = 1.2 * dx_left
        enddo

        ! right of origin
        do i=1, n_right
            ! position
            x(i + n_left) = 0.0 + dx_right*(i - 0.5)
            ! mass
            m(i + n_left) = rho_right*dx_right
            ! smoothing length
            h(i + n_left) = 1.2 * dx_right
        enddo

        ! mirror across -0.5 to make it periodic
        do i=1, n_left
            ! position
            x(n + i) = x_min - dx_left*(i - 0.5)
            ! mass
            m(n + i) = rho_left*dx_left
            ! smoothing length
            h(n + i) = 1.2 * dx_left
        enddo

        do i=1, n_right
            ! position
            x(n + i + n_left) = 0.0 + 2*x_min - dx_right*(i - 0.5)
            ! mass
            m(n + i + n_left) = rho_right*dx_right
            ! smoothing length
            h(n + i + n_left) = 1.2 * dx_right
        enddo

        ! setup velocity
        do i = 1, 2 * n
            v(i) = 0
        enddo

        ! after mirroing, we've doubled the number of particles, and changed xmin:
        n = n*2
        x_min = x_min - (x_max - x_min)

    end subroutine isothermal_setup


    subroutine sod_setup(x, v, rho, u, P, m, h, x_min, x_max, n_max, n, gamma)

        integer, intent(in) :: n_max
        integer :: n_left, n_right
        integer, intent(out) :: n
        real, parameter :: rho_left = 1.0, rho_right = 0.125, P_left = 1.0, P_right = 0.1
        real, intent(inout) :: x_min, x_max
        real :: dx_left = 0.001, dx_right = 0.008
        real, intent(in) :: gamma
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), P(n_max), u(n_max)

        ! calculate number of particles left and right of origin
        ! (round to nearest whole number)
        n_left = nint(abs(0 - x_min)/dx_left)
        n_right = nint(abs(x_max - 0)/dx_right)
        n = n_left + n_right

        ! left of origin
        do i=1, n_left
            ! position
            x(i) = x_min + dx_left*(i - 0.5)
            ! mass
            m(i) = rho_left*dx_left
            ! smoothing length
            h(i) = 1.2 * dx_left
            ! density
            rho(i) = rho_left
            ! Pressure
            P(i) = P_left
            ! internal energy according to adiabatic eos
            u(i) = P_left/((gamma - 1)*rho_left)
        enddo

        ! right of origin
        do i=1, n_right
            ! position
            x(i + n_left) = 0.0 + dx_right*(i - 0.5)
            ! mass
            m(i + n_left) = rho_right*dx_right
            ! smoothing length
            h(i + n_left) = 1.2 * dx_right
            ! density
            rho(i + n_left) = rho_right
            ! Pressure
            P(i + n_left) = P_right
            ! internal energy according to adiabatic eos
            u(i + n_left) = P_right/((gamma - 1)*rho_right)
        enddo

        ! mirror across -0.5 to make it periodic
        do i=1, n_left
            ! position
            x(n + i) = x_min - dx_left*(i - 0.5)
            ! mass
            m(n + i) = rho_left*dx_left
            ! smoothing length
            h(n + i) = 1.2 * dx_left
            ! density
            rho(n + i) = rho_left
            ! Pressure
            P(n + i) = P_left
            ! internal energy according to adiabatic eos
            u(n + i) = P_left/((gamma - 1)*rho_left)
        enddo

        do i=1, n_right
            ! position
            x(n + i + n_left) = 0.0 + 2*x_min - dx_right*(i - 0.5)
            ! mass
            m(n + i + n_left) = rho_right*dx_right
            ! smoothing length
            h(n + i + n_left) = 1.2 * dx_right
            ! density
            rho(n + i + n_left) = rho_right
            ! Pressure
            P(n + i + n_left) = P_right
            ! internal energy according to adiabatic eos
            u(n + i + n_left) = P_right/((gamma - 1)*rho_right)
        enddo

        ! setup velocity
        do i = 1, 2 * n
            v(i) = 0
        enddo

        ! after mirroing, we've doubled the number of particles, and changed xmin:
        n = n*2
        x_min = x_min - (x_max - x_min)
        
    end subroutine sod_setup

end module init
