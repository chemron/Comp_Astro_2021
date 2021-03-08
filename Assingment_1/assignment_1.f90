program assingment_1
    use derivs
    use evolution
    use init

    implicit none
    integer, parameter :: n_max = 150, n_ghosts = 10
    real, parameter :: x_min = 0.0, x_max = 1.0
    real, parameter :: c_0 = 1
    integer :: n
    real :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), &
    t_start, t_end, dt, dtout
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)
    print*, 'Hello World'

    ! initialise
    call setup(x, v, m, h, c_0, x_min, x_max, n_max, n)

    call set_ghosts(x, v, a, m, h, rho, u, P, c, x_min, x_max, n_max, n, n_ghosts)

    call get_derivs(x, v, a, m, h, rho, u, P, c, c_0, x_min, x_max, n_max, n, n_ghosts)

    call output(0.0, x, v, a, m, h, rho, u, P, c, n_max, n, n_ghosts, 0)

    ! evolve for 5 seconds
    t_start = 0.0
    t_end = 5.0
    dtout = 0.05
    dt = 0.1 * 0.2 * h(1) / c(1)

    call timestepping(x, v, a, m, h, rho, u, P, c, c_0, t_start, t_end, dt, dtout, x_min, x_max, n_max, n_ghosts, n)


end program assingment_1
