module evolution
    use init
    use derivs
    use outputs
    implicit none


contains

    subroutine leapfrog(x, v, a, m, h, rho, u, P, c, c_0, t, dt, x_min, x_max, n_max, n, n_ghosts)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: x_min, x_max, c_0, dt
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max), t
        real :: v_s(n_max), a_0(n_max)



        ! initialise
        a_0 = a

        x = x + dt * v + 0.5 * dt**2 * a_0
        v_s = v + dt * a_0
        call get_derivs(x, v_s, a, m, h, rho, u, P, c, c_0, x_min, x_max, n_max, n, n_ghosts)

        v = v_s + 0.5 * dt * (a - a_0)

        t = t + dt

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


    subroutine timestepping(x, v, a, m, h, rho, u, P, c, ke, c_0, t_start, t_end, dt, dtout, x_min, x_max, n_max, n_ghosts, n)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: x_min, x_max, c_0, t_start, t_end, dt, dtout
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), &
        h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max), ke(n_max)
        real :: t, tprint
        integer :: ifile

        t = t_start
        ifile = 1
        tprint = ifile * dtout
        do while (t <= t_end)
            call leapfrog(x, v, a, m, h, rho, u, P, c, c_0, t, dt, x_min, x_max, n_max, n, n_ghosts)
            ! calculate kinetic energy for every particle
            call get_kinetic_energy(v, m, ke, n, n_max)

            ! print the current total kinetic energy (sum over all particles)
            call print_ke(t, sum(ke))

            ! print output
            if (t >= tprint) then
                call output(t, x, v, a, m, h, rho, u, P, c, ke, n_max, n, n_ghosts, ifile)
                ifile = ifile + 1
                tprint = ifile * dtout
            endif
        enddo


    end subroutine timestepping




end module evolution