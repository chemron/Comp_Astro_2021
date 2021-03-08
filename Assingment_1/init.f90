module init
    implicit none
    real , parameter :: pi = 4.* atan (1.0)

    contains
    subroutine setup(x, v, m, h, c_0, x_min, x_max, n_max, n)

        integer, intent(in) :: n_max
        integer, intent(out) :: n
        real, parameter :: rho_0 = 1.0
        real, intent(in) :: x_min, x_max
        real :: dx
        real, intent(in) :: c_0
        integer :: i
        ! real, intent(in) :: x(n), v(n), m(n), h(n), rho(n), u(n), P(n), c(n)
        real, intent(out) :: x(n_max), v(n_max), m(n_max), h(n_max)
        n = 100

        ! setup position
        ! put particles in the middle of cells with width dx 
        dx = (x_max - x_min)/n
        do i=1, n
            x(i) = x_min + dx*(i - 0.5) 
        enddo

        ! setup velocity
        do i = 1, n
            v(i) = 10.0 ** (-4) * c_0 * sin(2 * pi * (x(i)-x_min) / (x_max-x_min))
        enddo

        ! setup mass
        do i = 1, n
            m(i) = rho_0*dx
        enddo

        ! setup smoothing length
        do i = 1, n
            h(i) = 1.2 * dx
        enddo


    end subroutine setup

    subroutine set_ghosts(x, v, a, m, h, rho, u, P, c, x_min, x_max, n_max, n, n_ghosts)

        integer, intent(in) :: n_max, n_ghosts, n
        real, intent(in) :: x_min, x_max
        integer :: i
        real, parameter :: rho_0 = 1.0
        real, intent(inout) :: x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max)
        real :: width

        ! get width
        width = abs(x_max - x_min)

        ! right ghosts
        do i=1, n_ghosts/2
            x(n + i) = x(i) + width 
            v(n + i) = v(i)
            a(n + i) = a(i)
            m(n + i) = m(i)
            h(n + i) = h(i)
            rho(n + i) = rho(i)
            u(n + i) = u(i)
            P(n + i) = P(i)
            c(n + i) = c(i)
        enddo
        
        ! left ghosts
        do i= 1, n_ghosts/2
            x(i + n + n_ghosts/2) = x(n + 1 - i) - width
            v(i + n + n_ghosts/2) = v(n + 1 - i)
            m(i + n + n_ghosts/2) = m(n + 1 - i)
            h(i + n + n_ghosts/2) = h(n + 1 - i)
            rho(i + n + n_ghosts/2) = rho(n + 1 - i)
            u(i + n + n_ghosts/2) = u(n + 1 - i)
            P(i + n + n_ghosts/2) = P(n + 1 - i)
            c(i + n + n_ghosts/2) = c(n + 1 - i)

        enddo

    end subroutine set_ghosts


    subroutine output(t, x, v, a, m, h, rho, u, P, c, n_max, n, n_ghosts, ifile)
        integer, intent(in) :: n_max, n, n_ghosts, ifile
        real, intent(in) :: t, x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max)
        integer :: lu = 1, i
        character(len=128) :: filename
    
        write(filename,"(a,i5.5)") 'output/snap_', ifile
    
        print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t    

        open(lu , file=filename, status='replace', action='write')
        write(lu,*) '# x, v, a, m, h, rho, u, P, c'
        write(lu,*) t
        do i=1,n + n_ghosts
            write(lu,*) x(i), v(i), a(i), m(i), h(i), rho(i), u(i), P(i), c(i)
        enddo
        close(lu)

    end subroutine output

end module init
