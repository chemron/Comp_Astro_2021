program sod
    use derivs
    use outputs
    use evolution
    use edges
    use init

    implicit none
    
    integer, parameter :: n_max = 1200, n_bound = 6
    real :: x_min = -0.5, x_max = 0.5, gamma = 1.4
    real, parameter :: c_0 = 1.0, alpha = 0.0, beta = 0.0
    integer :: n
    real :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), ke(n_max), &
    dudt(n_max), t_start, t_end, dt, dtout
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)

    ! initialise
    call sod_setup(x, v, rho, u, P, m, h, x_min, x_max, n_max, n, gamma)

    call get_derivs(x, v, a, m, h, rho, u, P, c, dudt, c_0, gamma, x_min, x_max, n_max, n, alpha, beta)
    
    call set_boundary(v, n, n_max, n_bound)

    call output(0.0, x, v, a, m, h, rho, u, P, c, ke, dudt, n_max, n, 0)

    ! initiallise ke file
    call initialise_ke_output()

    ! evolve for 0.1 seconds
    t_start = 0.0
    t_end = 0.2
    dtout = 0.002

    dt = 0.2 * minval(h(1:n)/c(1:n))
    print*, dt, x(n), v(n), a(n), m(n), h(n), rho(n), u(n), P(n), c(n)

    call timestepping(x, v, a, m, h, rho, u, P, c, dudt, ke, c_0, gamma, alpha, beta, t_start, t_end, dt, &
    dtout, x_min, x_max, n_max, n, n_bound)

end program sod
