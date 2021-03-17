module edges
    implicit none

    contains

    subroutine set_ghosts(x, v, m, h, rho, u, x_min, x_max, n_max, n, n_ghost)

        integer, intent(in) :: n_max, n
        integer, intent(out) :: n_ghost
        real, intent(in) :: x_min, x_max
        integer :: i
        real, intent(inout) :: x(n_max), v(n_max), m(n_max), h(n_max), rho(n_max), u(n_max)
        real :: width

        ! get width
        width = abs(x_max - x_min)

        n_ghost = 0

        do i = 1, n
            if (x(i) - 2 * h(i) .LE. x_min) then
                ! right ghosts
                n_ghost = n_ghost + 1
                x(n + n_ghost) = x(i) + width
                v(n + n_ghost) = v(i)
                m(n + n_ghost) = m(i)
                h(n + n_ghost) = h(i)
                rho(n + n_ghost) = rho(i)
                u(n + n_ghost) = u(i)
            elseif (x(i) + 2 * h(i) .GE. x_max) then
                ! left ghosts
                n_ghost = n_ghost + 1
                x(n + n_ghost) = x(i) - width
                v(n + n_ghost) = v(i)
                m(n + n_ghost) = m(i)
                h(n + n_ghost) = h(i)
                rho(n + n_ghost) = rho(i)
                u(n + n_ghost) = u(i)
            endif
        enddo

    end subroutine set_ghosts

    subroutine set_boundary(v, n, n_max, n_bound)
        integer, intent(in) :: n, n_max, n_bound
        real, intent(inout) :: v(n_max)
        integer :: i, n_mid

        do i = 1, n_bound
            ! left boundary:
            v(i) = 0.0
            ! right boundary:
            v(n - (i-1)) = 0.0
            
            n_mid = n/2
            ! middle left boundary:
            v(n_mid + i) = 0.0
            ! middle right boundary:
            v(n_mid - (i-1)) = 0.0
        enddo

    end subroutine set_boundary

end module edges