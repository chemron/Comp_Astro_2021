program main
    use output
    use integrators

    implicit none


    real :: x(2), v(2), t, t_print, dt, t_end, dt_out, ec, e, L
    integer :: i_file

    ! eccentricity
    ec = 0.7

    dt = 0.01
    ! 5000 timesteps, 1000 prints
    dt_out = dt * 5.0
    t_end = 5000.0 * dt
    t = 0.0

    ! Initial conditions:
    x(1) = 1.0 - ec
    x(2) = 0.0
    v(1) = 0.0
    v(2) = sqrt((1.0+ec)/(1.0-ec))

    i_file = 1
    t_print = i_file * dt_out

    ! make output file
    call initialise_output()

    ! Timestepping:
    do while (t <= t_end)

        ! single leapfrog step
        ! call leapfrog(x, v, dt)
        ! or, runge kutta:
        call rk4(x, v, dt)

        t = t + dt
        
        ! print output
        if (t >= t_print) then
            call get_conserved_quantities(x, v, L, e)
            call print_single_output(t, x, v, i_file)
            call print_output(t, x, v, e, L)
            i_file = i_file + 1
            t_print = i_file * dt_out
        endif

    enddo


end program main