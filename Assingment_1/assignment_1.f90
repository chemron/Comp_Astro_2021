program assingment_1
    use physics
    use init

    implicit none
    integer, parameter :: n_max = 150, n_ghosts = 10
    real, parameter :: c_0 = 1
    integer :: n
    real :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
    ! position (x), velocity (v), mass (m), smoothing length (h), density (rho), internal
    ! energy (u), pressure (P), and sound speed (c)
    print*, 'Hello World'


    call setup(x, v, m, h, c_0, n_max, n)

    call set_ghosts(x, v, m, h, n_max, n_ghosts, n)

    call get_density(x, m, h, rho, n_max, n_ghosts, n)

    call equation_of_state(rho, P, c, c_0, n_max, n, n_ghosts)

    call output(x, v, m, h, rho, u, P, c, n_max, n, n_ghosts)

end program assingment_1
