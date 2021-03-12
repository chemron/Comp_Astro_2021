module physics
    implicit none

contains

    subroutine get_kinetic_energy(v, m, ke, n, n_max)
        integer :: i
        integer, intent(in) :: n, n_max
        real, intent(in) :: v(n_max), m(n_max)
        real, intent(out) :: ke(n_max)

        do i = 1, n
            ke(i) = 0.5 * m(i) * v(i)**2
        enddo

    end subroutine get_kinetic_energy


    subroutine  get_smoothing_length(m, rho, h, n, n_max)
        integer :: i
        integer, parameter :: n_dim = 1
        real :: eta
        integer, intent(in) :: n, n_max
        real, intent(in) :: m(n_max), rho(n_max)
        real, intent(inout) :: h(n_max)
        ! TODO: justify eta/h_fac
        eta = 1.2
        do i = 1, n
            h(i) = eta * (m(i) / rho(i))**(1/n_dim)
        enddo

    end subroutine  get_smoothing_length


end module physics