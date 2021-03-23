module output
    implicit none

contains

    subroutine initialise_output()
        integer :: lu = 1

        open(lu , file='output.out', status='replace', action='write')
        write(lu,*) '# Time, x, y, dx, dy, e, L'
        close(lu)

    end subroutine initialise_output


    subroutine print_output(t, x, v, e, L)
        real, intent(in) :: t, x(2), v(2), e, L
        integer :: lu = 1
        open(lu , file='output.out', status='old', access='append')
        write(lu,*) t, x(1), x(2), v(1), v(2), e, L

        close(lu)

    end subroutine print_output


    subroutine print_single_output(t, x, dx, ifile)
        integer, intent(in) :: ifile
        real, intent(in) :: t, x(2), dx(2)
        integer :: lu = 1
        character(len=128) :: filename

        write(filename,"(a,i5.5)") 'output/snap_', ifile
        print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t
        open(lu , file=filename, status='replace', action='write')
        write(lu,*) '# x, y, dx, dy'
        write(lu,*) t
        write(lu,*) x(1), x(2), dx(1), dx(2)
        close(lu)

    end subroutine print_single_output

end module output
