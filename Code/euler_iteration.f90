
    subroutine euler_iteration(av,g)

!   This subroutine calculates the fluxes into each cell and then sums them to
!   update the primary flow properties

!   Explicitly declare the required variables
    use types
    use flux_stencil
    use smooth_stencil
    use flux_calculation
    implicit none
    type(t_appvars), intent(in) :: av
    type(t_grid), intent(inout) :: g
    real, dimension(g%ni,g%nj-1) :: mass_i, flux_i
    real, dimension(g%ni-1,g%nj) :: mass_j, flux_j
    integer :: i, j, ni, nj
    integer, parameter :: solution_order = 2

!   Get the block size and store locally for convenience
    ni = g%ni; nj = g%nj

!   Setup the continuity equation by calculating the mass flow through
!   the facets in both the i and j-directions. Store these values in
!   "mass_i" and "mass_j"
!   INSERT

    if (solution_order == 2) then
        call calculate_mass_flux_second_order(g, mass_i, mass_j)
    end if
    
!   END INSERT
     
!   Apply the wall boundary condition by checking that two nodes at the
!   end of a facet are both on a wall, if so then the appropriate mass
!   flow array is set to have zero flow through that facet
    where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
    where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!   Update the density with mass fluxes by calling "sum_fluxes"
!   INSERT

    call sum_fluxes(av, mass_i, mass_j, g%area, g%ro, g%dro)

!   END INSERT

!   Setup the conservation of energy equation by calculated the enthalpy flux
!   and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!   and "mass_j" from before
!   INSERT

    if (solution_order == 2) then
        call calculate_energy_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    end if

!   END INSERT

!   Update the internal energy with enthalpy fluxes
!   INSERT

    call sum_fluxes(av, flux_i, flux_j, g%area, g%roe, g%droe)

!   END INSERT

!   Setup the x-momentum equation including momentum flux and pressure forces
!   INSERT

    if (solution_order == 2) then
        call calculate_x_momentum_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    end if
    
!   END INSERT

!   Update the x-momentum with momentum flux
!   INSERT

    call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx, g%drovx)

!   END INSERT

!   Setup the y-momentum equation including momentum flux and pressure forces
!   INSERT

    if (solution_order == 2)then
        call calculate_y_momentum_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    end if

!   END INSERT

!   Update the y-momentum with momentum flux
!   INSERT

    call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy, g%drovy)

!   END INSERT

!   Add artificial viscosity by smoothing all of the primary flow variables
    call smooth_array(av,g%ro,g%corr_ro)
    call smooth_array(av,g%roe,g%corr_roe)
    call smooth_array(av,g%rovx,g%corr_rovx)
    call smooth_array(av,g%rovy,g%corr_rovy)
    

    end subroutine euler_iteration


