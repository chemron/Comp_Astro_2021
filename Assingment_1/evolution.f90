module evolution
    use init
    use derivs
    use outputs
    implicit none


contains

    subroutine leapfrog(x, v, a, m, h, rho, u, P, c, dudt, c_0, gamma, alpha, beta, t, dt, &
        x_min, x_max, n_max, n, n_ghosts, adiabatic)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: x_min, x_max, c_0, gamma, alpha, beta
        logical, intent(in) :: adiabatic
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max), dudt(n_max), t, dt
        real :: v_s(n_max), u_s(n_max), a_0(n_max), dudt_0(n_max)



        ! keep initial value
        a_0(1:n) = a(1:n)
        dudt_0(1:n) = dudt(1:n)

        ! update position
        x(1:n) = x(1:n) + dt * v(1:n) + 0.5 * dt**2 * a_0(1:n)
    
        ! get intermediate volicity
        v_s(1:n) = v(1:n) + dt * a_0(1:n)
        ! get inetmediat internal energy
        u_s(1:n) = u(1:n) + dt * dudt_0(1:n)

        call get_derivs(x, v, a, m, h, rho, u, P, c, dudt, c_0, gamma, x_min, x_max, n_max, &
        n, n_ghosts, adiabatic, alpha, beta)

        v(1:n) = v_s(1:n) + 0.5 * dt * (a(1:n) - a_0(1:n))
        u(1:n) = u_s(1:n) + 0.5 * dt * (dudt(1:n) - dudt_0(1:n))

        t = t + dt
        dt = 0.2 * minval(h(1:n)/c(1:n))


    end subroutine leapfrog


    subroutine get_kinetic_energy(v, m, ke, n, n_max)
        integer :: i
        integer, intent(in) :: n, n_max
        real, intent(in) :: v(n_max), m(n_max)
        real, intent(out) :: ke(n_max)

        do i = 1, n
            ke(i) = 0.5 * m(i) * v(i)**2
        enddo

    end subroutine get_kinetic_energy


    subroutine timestepping(x, v, a, m, h, rho, u, P, c, dudt, ke, c_0, gamma, alpha, beta, t_start, t_end, dt, &
        dtout, x_min, x_max, n_max, n_ghosts, n, n_bound, adiabatic)
        integer, intent(in) :: n, n_max, n_ghosts, n_bound
        real, intent(in) :: x_min, x_max, c_0, t_start, t_end, dtout, gamma, alpha, beta
        logical, intent(in) :: adiabatic
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), &
        h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), ke(n_max), dudt(n_max), dt
        real :: t, tprint
        integer :: ifile

        t = t_start
        ifile = 1
        tprint = ifile * dtout
        do while (t <= t_end)

            ! single leapfrog step
            call leapfrog(x, v, a, m, h, rho, u, P, c, dudt, c_0, gamma, alpha, beta, t, &
            dt, x_min, x_max, n_max, n, n_ghosts, adiabatic)

            ! set boundary conditions
            call set_boundary(v, n, n_max, n_bound)

            ! calculate kinetic energy for every particle
            call get_kinetic_energy(v, m, ke, n, n_max)

            ! print the current total kinetic energy (sum over all particles)
            call print_ke(t, sum(ke))

            ! print output
            if (t >= tprint) then
                call output(t, x, v, a, m, h, rho, u, P, c, ke, dudt, n_max, n, n_ghosts, ifile)
                ifile = ifile + 1
                tprint = ifile * dtout
            endif
        enddo


    end subroutine timestepping




end module evolution