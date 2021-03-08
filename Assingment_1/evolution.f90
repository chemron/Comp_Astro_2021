module evolution
    use init
    use derivs
    implicit none


contains

    subroutine leapfrog(x, v, a, m, h, rho, u, P, c, c_0, t, dt, x_min, x_max, n_max, n, n_ghosts)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: x_min, x_max, c_0, dt
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max), t
        real :: x_0(n_max), x_1(n_max), v_0(n_max), v_half(n_max), v_1(n_max), a_0(n_max), a_1(n_max)  

        ! initialise
        x_0 = x
        v_0 = v
        a_0 = a


        x_1 = x_0 + dt * v_0 + 0.5 * dt**2 * a_0
        v_half = v_0 + dt * a_0
        call get_derivs(x, v, a, m, h, rho, u, P, c, c_0, x_min, x_max, n_max, n, n_ghosts)

        v_1 = v_half + 0.5 * dt * (a_1 - a_0)

        ! finalise

        x = x_1
        v = v_1
        a = a_1
        t = t + dt

    end subroutine leapfrog


    subroutine timestepping(x, v, a, m, h, rho, u, P, c, c_0, t_start, t_end, dt, dtout, x_min, x_max, n_max, n_ghosts, n)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: x_min, x_max, c_0, t_start, t_end, dt, dtout
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), u(n_max), P(n_max), c(n_max)
        real :: t, tprint
        integer :: ifile
        
        t = t_start
        ifile = 1
        tprint = ifile * dtout
        do while (t <= t_end)
            call leapfrog(x, v, a, m, h, rho, u, P, c, c_0, t, dt, x_min, x_max, n_max, n, n_ghosts)
            if (t >= tprint) then
                call output(t, x, v, a, m, h, rho, u, P, c, n_max, n, n_ghosts, ifile)
                ifile = ifile + 1
                tprint = ifile * dtout
            endif
        enddo

    end subroutine timestepping




end module evolution