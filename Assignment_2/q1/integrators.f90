module integrators
 
contains
    ! we have: 
    ! x1  = x
    ! x1' = f1(x2) = x2  ( = v )
    ! x2' = f2(x1) = -(x1)/(norm2(x1)**3) ( = a )
    

    ! velocity
    function f1(x2)
        real, intent(in) :: x2(2)
        real :: f1(2)
        ! I know it's super obvious i.e. velocity = velocity, but just wanted it
        ! for completeness
        f1 = x2
    end function f1


    ! acceleration
    function f2(x1)
        ! position
        real, intent(in) :: x1(2)
        real :: f2(2), r

        r = norm2(x1)
        f2 = -x1/(r**3)

    end function f2


    subroutine leapfrog(x, v, dt)
        real, intent(inout) :: x(2), v(2)
        real, intent(in) :: dt
        real :: a_1(2), a_2(2), v_s(2)
  
        ! get initial acceleration
        a_1 = f2(x)

        ! update position
        x = x + dt*v + 0.5*dt**2*a_1

        ! get intermediate velocity
        v_s = v + dt * a_1

        ! update acceleration
        a_2 = f2(x)

        ! update velocity
        v = v_s + 0.5 * dt * (a_2 - a_1)

    end subroutine leapfrog


    ! 4th order runge-kutta method
    ! x1 = x, x2 = v
    subroutine rk4(x1, x2, dt)
        real, intent(inout) :: x1(2), x2(2)
        real, intent(in) :: dt
        real :: k11(2), k12(2), k13(2), k14(2), k21(2), k22(2), k23(2), k24(2)

        k11 = dt * f1(x2)
        k21 = dt * f2(x1)

        k12 = dt * f1(x2 + 0.5 * k21)
        k22 = dt * f2(x1 + 0.5 * k11)

        k13 = dt * f1(x2 + 0.5 * k22)
        k23 = dt * f2(x1 + 0.5 * k12)

        k14 = dt * f1(x2 + k23)
        k24 = dt * f2(x1 + k13)

        ! update position
        x1 = x1 + (1.0/6.0) * (k11 + 2*k12 + 2*k13 + k14) 
        ! update velocity
        x2 = x2 + (1.0/6.0) * (k21 + 2*k22 + 2*k23 + k24)

    end subroutine rk4


    subroutine get_conserved_quantities(x, v, L, e)
        real, intent(in) :: x(2), v(2)
        real, intent(out) :: L, e
        real :: r 

        r = norm2(x)
        L = x(1)*v(2) - x(2)*v(1)
        e = 0.5*(v(1)**2 + v(2)**2) - 1/r

    end subroutine get_conserved_quantities


end module integrators

