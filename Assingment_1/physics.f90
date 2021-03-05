module physics
    implicit none

contains
    real function W(x_a, x_b, h)
        real, intent(in) :: x_a, x_b, h
        real :: q
        integer, parameter :: d=1

        q = abs(x_a - x_b)/h

        if ((0.0 <= q) .and. (q < 1)) then
            w = (1.0/4.0)*(2.0-q)**3 - (1 - q)**3
        elseif ((1.0 <= q) .and. (q < 2.0)) then
            w = (1.0/4.0)*(2.0-q)**3
        else
            w = 0
        endif

        W = (1.0/h**d) * (2.0/3.0) * w


        
    end function W

    subroutine get_density(x, m, h, rho, n_max, n_ghosts, n)
        integer, intent(in) :: n_max, n_ghosts, n
        integer :: a, b
        real, intent(in) :: x(n_max), m(n_max), h(n_max)
        real, intent(out) :: rho(n_max)

        do a = 1, n + n_ghosts
            rho(a) = 0
            ! summation:
            do b = 1, n + n_ghosts
                rho(a) = rho(a) + m(b) * W(x(a), x(b), h(a))
            enddo

        enddo

    end subroutine get_density

end module physics
