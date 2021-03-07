module evolution
    implicit none


contains

    subroutine leapfrog(x, v, a, m, h, rho, P, c, c_0, dt, n_max, n_ghosts, n)
        integer, intent(in) :: n, n_max, n_ghosts
        real, intent(in) :: m(n_max), h(n_max), rho(n_max), P(n_max), c(n_max), c_0, dt
        real, intent(inout) :: x(n_max), v(n_max), a(n_max)
        real :: x_0(n_max), x_1(n_max), v_0(n_max), v_half(n_max), v_1(n_max), a_0(n_max), a_1(n_max)  

        ! initialise
        x_0 = x
        v_0 = v
        a_0 = a


        x_1 = x_0 + dt * v_0 + 0.5 * dt**2 * a_0
        v_half = v_0 + dt * a_0
        call get_derivs(x_1, a_1, m, h, rho, P, c, c_0, n_max, n_ghosts, n)

        v_1 = v_half + 0.5 * dt * (a_1 - a_0)

        ! finalise

        x = x_1
        v = v_1
        a = a_1

    end subroutine leapfrog

end module evolution