subroutine flow_guess(av,g,bcs,guesstype)

!   Initialise conserved variables using the selected boundary mode and guess.
    use types
    implicit none
    type(t_appvars), intent(inout) :: av
    type(t_grid), intent(inout) :: g
    type(t_bconds), intent(in) :: bcs
    integer, intent(in) :: guesstype
    integer :: i, j, ni, nj, j_mid
    real :: ro_in, t_in, v_in, lx, ly, l
    real :: t_out, v_out, ro_out, m_dot_out
    real :: m_dot, mach_lim, t_lim
    real :: l_i(g%ni)
    real :: v_guess(g%ni), ro_guess(g%ni), t_guess(g%ni), p_guess(g%ni)

!   parameters for guesstype == 2
    real, parameter :: tol = 1e-6
    integer, parameter :: max_iter = 100
    integer :: iter
    real :: mach_lower, mach_upper, residual
    real :: nd_m_dot(g%ni)
    real :: mach_guess(g%ni), compressible_factor(g%ni)

    ni = g%ni; nj = g%nj

    if (av%supersonic == 0) then
!       Restore master's outlet-based subsonic initialization.
        t_out = bcs%tstag * (bcs%p / bcs%pstag)**av%fgam
        v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
        ro_out = bcs%p / (av%rgas * t_out)

        if (guesstype == 1) then
            g%ro = ro_out
            g%roe = g%ro * (av%cv * t_out + 0.5 * v_out**2)
            j_mid = nj / 2
            do i = 1, ni-1
                lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid)
                l = hypot(lx,ly)
                g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
                g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
            end do
            g%rovx(ni,:) = g%rovx(ni-1,:)
            g%rovy(ni,:) = g%rovy(ni-1,:)
            write(6,*) 'Crude subsonic flow guess calculated'
        else if (guesstype == 2) then
            l_i = hypot(g%x(:,nj) - g%x(:,1), g%y(:,nj) - g%y(:,1))
            m_dot_out = ro_out * v_out * l_i(ni)
            mach_lim = 0.99
            t_lim = bcs%tstag * (1.0 + ((av%gam - 1.0) / 2.0) * mach_lim**2)**(-1.0)
            v_guess = m_dot_out/(ro_out*l_i)
            t_guess = bcs%tstag - v_guess**2 / (2*av%cp)
            t_guess = max(t_guess, t_lim)
            p_guess = bcs%pstag * (t_guess/bcs%tstag)**(1.0/av%fgam)
            ro_guess = p_guess/(av%rgas * t_guess)
            v_guess = m_dot_out/(ro_guess*l_i)
            do j = 1, nj
                do i = 1, ni-1
                    lx = g%lx_j(i,j); ly = g%ly_j(i,j)
                    l = hypot(lx,ly)
                    g%ro(i,j) = ro_guess(i)
                    g%roe(i,j) = ro_guess(i) * (av%cv * t_guess(i) + 0.5*v_guess(i)**2)
                    g%rovx(i,j) = ro_guess(i) * v_guess(i) * ly / l
                    g%rovy(i,j) = -ro_guess(i) * v_guess(i) * lx / l
                end do
            end do
            g%ro(ni,:) = g%ro(ni-1,:)
            g%roe(ni,:) = g%roe(ni-1,:)
            g%rovx(ni,:) = g%rovx(ni-1,:)
            g%rovy(ni,:) = g%rovy(ni-1,:)
            write(6,*) 'Refined subsonic flow guess calculated'
        else
            error stop 'Unknown flow guesstype; expected 1 or 2'
        end if
    else

    !   Both supersonic guesses use the prescribed inlet state.
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



    !       generate a non-dimentional mass flowrate list
            l_i = hypot(g%x(:,nj) - g%x(:,1), g%y(:,nj) - g%y(:,1))
            m_dot = ro_in * v_in * l_i(1)
            nd_m_dot = m_dot * sqrt(av%cp * bcs%tstag)/(l_i * bcs%pstag)

    !       initialize a list of guess mach number using inlet mach number
            mach_guess = sqrt(2/(av%gam - 1) * (bcs%tstag/t_in - 1))

    !       Bisect the supersonic branch, where mass flow decreases with Mach.
            if (any(nd_m_dot <= 0.0) .or. any(nd_m_dot > calc_nd_m_dot(1.0) + tol)) then
                error stop 'Mass flow target has no supersonic Mach solution'
            end if
            do i = 1, ni
                mach_lower = 1.0
                mach_upper = max(2.0, mach_guess(i))
                do iter = 1, max_iter
                    if (calc_nd_m_dot(mach_upper) <= nd_m_dot(i)) exit
                    mach_upper = 2.0 * mach_upper
                end do
                if (iter > max_iter) error stop 'Could not bracket Mach solution'

                do iter = 1, max_iter
                    mach_guess(i) = 0.5 * (mach_lower + mach_upper)
                    residual = calc_nd_m_dot(mach_guess(i)) - nd_m_dot(i)
                    if (abs(residual) <= tol) exit
                    if (residual > 0.0) then
                        mach_lower = mach_guess(i)
                    else
                        mach_upper = mach_guess(i)
                    end if
                end do
                if (iter > max_iter) error stop 'Mach iteration did not converge'
            end do

    !       print mach number for inspection
            write(6,*) "inlet mach number: ", mach_guess(1), ", outlet mach number: ", mach_guess(ni)

            compressible_factor = 1+(av%gam - 1)/2 * mach_guess**2
            t_guess = bcs%tstag/compressible_factor
            p_guess = bcs%pstag/compressible_factor**(av%gam/(av%gam - 1))
            ro_guess = bcs%rostag/compressible_factor**(1/(av%gam-1))
            v_guess = m_dot/(ro_guess * l_i)

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

    contains

    elemental function calc_nd_m_dot(M)
        real, intent(in) :: M
        real :: calc_nd_m_dot

        calc_nd_m_dot = av%gam/sqrt(av%gam - 1) * M * &
            (1 + (av%gam - 1)/2 * M**2)**(-0.5*(av%gam + 1)/(av%gam - 1))
    end function calc_nd_m_dot

end subroutine flow_guess
