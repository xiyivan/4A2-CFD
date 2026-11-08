subroutine flow_guess(av,g,bcs,guesstype)

!   Initialise the conserved variables from the prescribed inlet state.
    use types
    implicit none
    type(t_appvars), intent(inout) :: av
    type(t_grid), intent(inout) :: g
    type(t_bconds), intent(in) :: bcs
    integer, intent(in) :: guesstype
    integer :: i, j, ni, nj
    real :: ro_in, t_in, v_in, lx, ly, l
    real :: m_dot, mach_lim, t_lim
    real :: l_i(g%ni)
    real :: v_guess(g%ni), ro_guess(g%ni), t_guess(g%ni), p_guess(g%ni)

    ni = g%ni; nj = g%nj

!   Both guesses use the same prescribed supersonic inlet state.
    t_in = bcs%tstag * (bcs%p / bcs%pstag)**av%fgam
    ro_in = bcs%p / (av%rgas * t_in)
    v_in = sqrt(2.0 * av%cp * (bcs%tstag - t_in))

    if (guesstype == 1) then

!       Uniform inlet state throughout the domain, including flow direction.
        g%ro = ro_in
        g%rovx = ro_in * v_in * cos(bcs%alpha)
        g%rovy = ro_in * v_in * sin(bcs%alpha)
        g%roe = ro_in * (av%cv * t_in + 0.5 * v_in**2)
        write(6,*) 'Uniform inlet flow guess calculated'

    else if (guesstype == 2) then

!       varying-flow guess based on continuity.
        l_i = hypot(g%x(:,nj) - g%x(:,1), g%y(:,nj) - g%y(:,1))
        m_dot = ro_in * v_in * l_i(1)

        mach_lim = 1.01
        t_lim = bcs%tstag / (1.0 + 0.5 * (av%gam - 1.0) * mach_lim**2)

        v_guess = m_dot / (ro_in * l_i)
        t_guess = min(t_lim, bcs%tstag - v_guess**2 / (2.0 * av%cp))
        if (any(t_guess <= 0.0)) then
            error stop 'Varying flow guess gives nonpositive temperature; use guesstype = 1'
        end if
        p_guess = bcs%pstag * (t_guess / bcs%tstag)**(1.0 / av%fgam)
        ro_guess = p_guess / (av%rgas * t_guess)
        v_guess = m_dot / (ro_guess * l_i)

!       Direct velocity along the j = const gridlines.
        do j = 1, nj
            do i = 1, ni-1
                lx = g%lx_j(i,j); ly = g%ly_j(i,j)
                l = hypot(lx,ly)
                g%ro(i,j) = ro_guess(i)
                g%roe(i,j) = ro_guess(i) * (av%cv * t_guess(i) + 0.5 * v_guess(i)**2)
                g%rovx(i,j) = ro_guess(i) * v_guess(i) * ly / l
                g%rovy(i,j) = -ro_guess(i) * v_guess(i) * lx / l
            end do
        end do

        g%ro(ni,:) = g%ro(ni-1,:)
        g%roe(ni,:) = g%roe(ni-1,:)
        g%rovx(ni,:) = g%rovx(ni-1,:)
        g%rovy(ni,:) = g%rovy(ni-1,:)
        write(6,*) 'Refined flow guess calculated'

    else
        error stop 'Unknown flow guesstype; expected 1 or 2'
    end if

    write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
        g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
    write(6,*)

!   Reference values used to normalise convergence measures.
    av%ro_ref = sum(g%ro(1,:)) / nj
    av%roe_ref = sum(g%roe(1,:)) / nj
    av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj


!   initialize corrective primary properties to zero
    g%corr_ro = 0
    g%corr_roe = 0
    g%corr_rovx = 0
    g%corr_rovy = 0
end subroutine flow_guess
