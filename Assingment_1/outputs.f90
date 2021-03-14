module outputs
    implicit none

contains

    subroutine initialise_ke_output()
        integer :: lu = 1

        open(lu , file='energy.out', status='replace', action='write')
        write(lu,*) '# t, ke'
        close(lu)

    end subroutine initialise_ke_output


    subroutine print_ke(t, ke)
        real, intent(in) :: t, ke
        integer :: lu = 1
        open(lu , file='energy.out', status='old', access='append')
        ! '# t, ke'
        write(lu,*) t, ke

        close(lu)

    end subroutine print_ke


    subroutine output(t, x, v, a, m, h, rho, u, P, c, ke, n_max, n, n_ghosts, ifile)
        integer, intent(in) :: n_max, n, n_ghosts, ifile
        real, intent(in) :: t, x(n_max), v(n_max), a(n_max), m(n_max), h(n_max), rho(n_max), &
        u(n_max), P(n_max), c(n_max), ke(n_max)
        integer :: lu = 1, i
        character(len=128) :: filename

        write(filename,"(a,i5.5)") 'output/snap_', ifile
        print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t
        open(lu , file=filename, status='replace', action='write')
        write(lu,*) '# x, v, a, m, h, rho, u, P, c, ke'
        write(lu,*) t
        do i=1,n + n_ghosts
            write(lu,*) x(i), v(i), a(i), m(i), h(i), rho(i), u(i), P(i), c(i), ke(i)
        enddo
        close(lu)

    end subroutine output

end module outputs
