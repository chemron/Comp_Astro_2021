program sod
    use derivs
    use outputs
    use evolution
    use init

    implicit none
    
    integer, parameter :: n_max = 1200, n_ghosts = 0, n_bound = 10
    real, parameter :: x_min = 0.0, x_max = 1.0, gamma = 1.4
    real, parameter :: c_0 = 1.0, alpha = 1.0, beta = 2.0
    logical, parameter :: adiabatic = .True.
    integer :: n
    real :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), ke(n_max), &
    dudt(n_max), t_start, t_end, dt, dtout
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)

    ! initialise
    call sod_setup(x, v, rho, u, P, m, h, n_max, n, gamma)

    call get_derivs(x, v, a, m, h, rho, u, P, c, dudt, c_0, gamma, x_min, x_max, n_max, n, n_ghosts, adiabatic, alpha, beta)
    
    call set_boundary(v, n, n_max, n_bound)

    call output(0.0, x, v, a, m, h, rho, u, P, c, ke, dudt, n_max, n, n_ghosts, 0)

    ! initiallise ke file
    call initialise_ke_output()

    ! evolve for 0.1 seconds
    t_start = 0.0
    t_end = 0.2
    dtout = 0.002
    dt = 0.2 * minval(h(1:n)/c(1:n))

    call timestepping(x, v, a, m, h, rho, u, P, c, dudt, ke, c_0, gamma, t_start, t_end, dt, &
    dtout, x_min, x_max, n_max, n_ghosts, n, n_bound, adiabatic)

end program sod
