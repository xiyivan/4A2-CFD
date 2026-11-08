    
    subroutine apply_bconds(av,g,bcs)

!   This subroutine applies both the inlet and outlet boundary conditions, as
!   it modifies both the primary and secondary flow variables they must be
!   calculated first

!   Explicitly declare the required variables
    use types
    implicit none
    type(t_appvars), intent(in) :: av
    type(t_grid), intent(inout) :: g
    type(t_bconds), intent(inout) :: bcs

!   Declare the other variables you need here
!   INSERT
    real :: v(g%nj), T(g%nj)

!   Calculate Mach number at the inlet
    T = bcs%tstag * (bcs%p/bcs%pstag)**(av%fgam)
    bcs%ro = bcs%p/(av%rgas * T)

!   Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!   "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!   "hstag(1,:)"
!   INSERT

    g%ro(1, :) = bcs%ro
    g%p(1, :) = bcs%p
    v = sqrt((bcs%tstag - T)*2*av%cp)
    g%vx(1,:) = v * cos(bcs%alpha)
    g%vy(1,:) = v * sin(bcs%alpha)
    g%rovx(1,:) = g%vx(1,:) * bcs%ro
    g%rovy(1,:) = g%vy(1,:) * bcs%ro
    g%roe(1,:) = bcs%ro * (av%cv * T + 0.5 * v**2)
    g%hstag(1,:) = av%cp * bcs%tstag

!   END INSERT


    end subroutine apply_bconds


