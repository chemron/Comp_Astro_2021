program isothermal
    use derivs
    use outputs
    use evolution
    use init

    implicit none
    
    integer, parameter :: n_max = 600, n_ghosts = 0
    real, parameter :: x_min = 0.0, x_max = 1.0
    real, parameter :: c_0 = 1
    integer :: n
    real :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), ke(n_max), &
    t_start, t_end, dt, dtout
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)

    ! initialise
    ! call setup(x, v, m, h, c_0, x_min, x_max, n_max, n)
    call isothermal_setup(x, v, m, h, c_0, n_max, n)

    call get_density(x, m, h, rho, n_max, n_ghosts, n)

    call iso_output(0.0, x, v, a, m, h, rho, u, P, c, ke, n_max, n, n_ghosts, 0)

end program isothermal
