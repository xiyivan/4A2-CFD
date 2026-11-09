    
    subroutine read_settings(av,bcs)

!   Read in the application variables and gas constants, they are in the
!   standard input file which is already assigned to unit 5

!   Explicitly declare the required variables
    use types
    implicit none
    type(t_appvars), intent(out) :: av
    type(t_bconds), intent(out) :: bcs
    character(len=80) :: tempname
    character(len=512) :: options
    character(len=80) :: extra
    integer :: ios


!   Read the case name and trim to the required length
    read(5,*) tempname
    av%casename = trim(tempname)

!   You should read in the following variables sequentially and store them in
!   the dervived "av" datatype with the % syntax:
!       gam, rgas 
!       cfl, sfac, d_max
!       nsteps
!       ni, nj
!   INSERT

    read(5,*) av%gam, av%rgas
    read(5, *) av%cfl, av%sfac, av%d_max
    read(5, *) av%nsteps
    read(5, *) av%ni, av%nj

 !    END INSERT

!   Calculate other gas constants used throughout the calculation
    av%cp = av%rgas * av%gam / (av%gam - 1.0)
    av%cv = av%cp / av%gam
    av%fgam = (av%gam - 1.0) / av%gam

!   Scale the smoothing factor and convergence limits by cfl number, the 
!   changes over a timestep should be proportional to dt, which is 
!   proportional to cfl
    av%sfac = av%sfac * av%cfl
    av%d_max = av%d_max * av%cfl

!   Average convergence limit on residuals is set to half of the maximum
    av%d_avg = 0.5 * av%d_max

!   Read the inlet boundary condition data and store into the "bcs" datatype
!       pstag, tstag, alpha, rfin
!   INSERT

    read(5,*) bcs%pstag, bcs%tstag, bcs%alpha, bcs%rfin

!   END INSERT

!   Convert the inlet angle to radians
    bcs%alpha = bcs%alpha * 3.14159 / 180.0

!   Calculate the inlet stagnation density "rostag"
!   INSERT

    bcs%rostag = bcs%pstag / (av%rgas * bcs%tstag)

!   END INSERT

!   Read static pressure (outlet for subsonic, inlet for supersonic mode)
!   INSERT

    read(5, *) bcs%p

!   An optional final line selects boundary mode, guess and correction strength.
!   EOF preserves defaults for legacy seven-line input files.
    read(5,'(A)',iostat=ios) options
    if (ios > 0) error stop 'Could not read optional feature settings'
    if (ios == 0) then
        read(options,*,iostat=ios) av%supersonic, av%guesstype, av%fcorr
        if (ios /= 0) error stop 'Expected feature settings: supersonic guesstype fcorr'
        read(options,*,iostat=ios) av%supersonic, av%guesstype, av%fcorr, extra
        if (ios == 0) error stop 'Expected only three feature settings'
    end if
    if (av%supersonic /= 0 .and. av%supersonic /= 1) error stop 'supersonic must be 0 or 1'
    if (av%guesstype /= 1 .and. av%guesstype /= 2) error stop 'guesstype must be 1 or 2'
    if (.not. (av%fcorr >= 0.0 .and. av%fcorr <= 1.0)) error stop 'fcorr must be between 0 and 1'

!   END INSERT

!   Print the settings to check they have been read, you can use this syntax
!   anywhere else you want in the program to debug your code
    write(6,*)
    write(6,*) 'Solver begins on ', av%casename, ' case'
    write(6,*)
    write(6,*) 'Read application variables from file'
    write(6,*) '  rgas =', av%rgas, 'cp =', av%cp, 'cv =', av%cv
    write(6,*) '  CFL =', av%cfl, 'sfac =', av%sfac
    write(6,*) '  Convergence  d_max =', av%d_max
    write(6,*) '  Features: supersonic =', av%supersonic, 'guesstype =', av%guesstype, 'fcorr =', av%fcorr
    if (av%supersonic == 1) then
        write(6,*) '  Static pressure specifies the supersonic inlet'
    else
        write(6,*) '  Static pressure specifies the subsonic outlet'
    end if
    write(6,*) '  Mesh size  ni =', av%ni, 'nj =', av%nj
    write(6,*) '  Inlet  pstag =', bcs%pstag, 'tstag =', bcs%tstag, &
        'alpha = ', bcs%alpha, 'p = ', bcs%p
    write(6,*)

    end subroutine read_settings


