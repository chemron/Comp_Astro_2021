program assingment_1
    use derivs
    use init

    implicit none
    integer, parameter :: n_max = 150, n_ghosts = 10
    real, parameter :: x_min = 0.0, x_max = 1.0
    real, parameter :: c_0 = 1
    integer :: n
    real :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)
    print*, 'Hello World'


    call setup(x, v, m, h, c_0, x_min, x_max, n_max, n)

    call set_ghosts(x, v, a, m, h, rho, u, P, c, x_min, x_max, n_max, n, n_ghosts)

    call get_derivs(x, v, a, m, h, rho, u, P, c, c_0, x_min, x_max, n_max, n, n_ghosts)

    call output(x, v, a, m, h, rho, u, P, c, n_max, n, n_ghosts)


end program assingment_1
