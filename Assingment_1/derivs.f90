module derivs
    use init
    use physics
    implicit none

contains

    subroutine kernal(x_a, x_b, h, W, grad_W)
        real, intent(in) :: x_a, x_b, h
        real, intent(out) :: W, grad_W
        real :: q, grad_q, dw
        integer, parameter :: d=1

        ! calculate q
        q = abs(x_a - x_b)/h

        ! calculate grad_q
        grad_q = (1.0/h) * (x_a - x_b)/ abs(x_a - x_b)

        ! calculate W

        if ((0.0 <= q) .and. (q < 1)) then
            w = (1.0/4.0)*(2.0-q)**3 - (1 - q)**3
        elseif ((1.0 <= q) .and. (q < 2.0)) then
            w = (1.0/4.0)*(2.0-q)**3
        else
            w = 0
        endif

        W = (1.0/h**d) * (2.0/3.0) * w

        ! calculate dw/dq

        if ((0.0 <= q) .and. (q < 1)) then
            dw = -(3.0/4.0)*(2.0-q)**2 + 3.0 * (1 - q)**2
        elseif ((1.0 <= q) .and. (q < 2.0)) then
            dw = -(3.0/4.0)*(2.0-q)**2
        else
            dw = 0
        endif

        ! grad_W = dw/dq * grad_q

        dW = (1.0/h**d) * (2.0/3.0) * dw

        grad_W = dW * grad_q

    end subroutine kernal


    subroutine get_density(x, m, h, rho, n_max, n_ghosts, n)
        integer, intent(in) :: n_max, n_ghosts, n
        integer :: a, b
        real :: W, grad_W
        real, intent(in) :: x(n_max), m(n_max), h(n_max)
        real, intent(out) :: rho(n_max)

        do a = 1, n
            rho(a) = 0
            ! summation:
            do b = 1, n + n_ghosts
                call kernal(x(a), x(b), h(a), W, grad_W)
                rho(a) = rho(a) + m(b) * W
            enddo

        enddo

    end subroutine get_density


    subroutine equation_of_state(rho, P, c, c_0, n_max, n, n_ghosts)
        integer, intent(in) :: n_max, n_ghosts, n
        integer :: i
        real, intent(in) :: rho(n_max), c_0
        real, intent(out) :: P(n_max), c(n_max)

        ! use eos from lecture 1
        do i = 1, n + n_ghosts
            c(i) = c_0
            P(i) = c(i)**2 * rho(i)
        enddo

    end subroutine equation_of_state


    subroutine viscosity(r_a, r_b, v_a, v_b, rho_a, rho_b, c_a, c_b, q_a, q_b, alpha, beta)
        real, intent(in) :: r_a, r_b, v_a, v_b, rho_a, rho_b, c_a, c_b, alpha, beta
        real, intent(out) :: q_a, q_b
        real :: r_ab, v_ab, v_sig_a, v_sig_b
        
        r_ab = (r_a - r_b)/abs(r_a - r_b)
        v_ab = (v_a - v_b)
        v_sig_a = alpha * c_a - beta * (v_ab * r_ab)
        v_sig_b = alpha * c_b - beta * (v_ab * r_ab)

        q_a = 0.0
        q_b = 0.0
        if (v_ab * r_ab < 0.0) then
            q_a = (-0.5) * rho_a * v_sig_a * v_ab * r_ab
            q_b = (-0.5) * rho_b * v_sig_b * v_ab * r_ab
        endif

    end subroutine viscosity


    subroutine get_accel(x, v, a, m, h, rho, P, c, n_max, n_ghosts, n, alpha, beta)
        integer, intent(in) :: n_max, n_ghosts, n
        integer :: i, j
        real :: W, grad_W_i, grad_W_j, q_i, q_j
        real, intent(in) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), P(n_max), c(n_max), alpha, beta
        real, intent(out) :: a(n_max)

        do i = 1, n
            a(i) = 0.0

            do j = 1, n + n_ghosts
                ! grad_W should be zero if i == j
                if (i /= j) then
                    call viscosity(x(i), x(j), v(i), v(j), rho(i), rho(j), c(i), c(j), q_i, q_j, alpha, beta)
                    call kernal(x(i), x(j), h(i), W, grad_W_i)
                    call kernal(x(i), x(j), h(j), W, grad_W_j)
                    a(i) = a(i) - m(j) * ((P(i) + q_i)/rho(i)**2 * grad_W_i + (P(j) + q_j)/rho(j)**2 * grad_W_j)
                endif
            enddo
        enddo

    end subroutine get_accel


    subroutine get_derivs(x, v, a, m, h, rho, u, P, c, c_0, x_min, x_max, n_max, n, n_ghosts)
        integer :: i
        integer, intent(in) :: n_max, n_ghosts, n
        real, intent(in) :: c_0, x_min, x_max
        real, parameter :: rho_0 = 1.0
        real :: alpha, beta
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max)
        alpha = 1
        beta = 2

        call set_ghosts(x, v, a, m, h, rho, u, P, c, x_min, x_max, n_max, n, n_ghosts)

        do i = 1, 3
            call get_density(x, m, h, rho, n_max, n_ghosts, n)
            call get_smoothing_length(m, rho, h, n, n_max)
        enddo
    

        call equation_of_state(rho, P, c, c_0, n_max, n, n_ghosts)


        ! TODO: remove pressure, soundspeed, acceleration
        call set_ghosts(x, v, a, m, h, rho, u, P, c, x_min, x_max, n_max, n, n_ghosts)
        

        call get_accel(x, v, a, m, h, rho, P, c, n_max, n_ghosts, n, alpha, beta)


    end subroutine get_derivs

end module derivs
