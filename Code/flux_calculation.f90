module flux_calculation
    use types
    implicit none

    
contains

subroutine calculate_mass_flux_second_order(g,  mass_i, mass_j)
    type(t_grid), intent(inout) :: g
    real, intent(out) ::  mass_i(:,:), mass_j(:,:)
    integer :: ni, nj
    
    ni = g%ni
    nj = g%nj

    mass_i = g%lx_i * 0.5 * (g%rovx(:, 1:nj-1) + g%rovx(:, 2:nj)) + &
            g%ly_i * 0.5 * (g%rovy(:, 1:nj-1) + g%rovy(:, 2:nj))
    mass_j = g%lx_j * 0.5 * (g%rovx(1:ni-1, :) + g%rovx(2:ni, :)) + &
            g%ly_j * 0.5 * (g%rovy(1:ni-1, :) + g%rovy(2:ni, :))
    
end subroutine calculate_mass_flux_second_order

subroutine calculate_mass_flux_fourth_order(g,  mass_i, mass_j)
    type(t_grid), intent(inout) :: g
    real, intent(out) ::  mass_i(:,:), mass_j(:,:)
    integer :: ni, nj
    
    ni = g%ni
    nj = g%nj

    
    
end subroutine calculate_mass_flux_fourth_order

subroutine calculate_energy_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    type(t_grid), intent(inout) :: g
    real, intent(in) :: mass_i(:,:), mass_j(:,:)
    real, intent(inout) :: flux_i(:,:), flux_j(:,:)
    integer :: ni, nj

    ni = g%ni
    nj = g%nj

    flux_i = mass_i * 0.5 * (g%hstag(:, 1:nj-1) + g%hstag(:, 2:nj))
    flux_j = mass_j * 0.5 * (g%hstag(1:ni-1, :) + g%hstag(2:ni, :))
    
end subroutine calculate_energy_flux_second_order

subroutine calculate_x_momentum_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    type(t_grid), intent(inout) :: g
    real, intent(in) :: mass_i(:,:), mass_j(:,:)
    real, intent(inout) :: flux_i(:,:), flux_j(:,:)
    integer :: ni, nj

    ni = g%ni
    nj = g%nj

    flux_i = mass_i * 0.5 * (g%vx(:, 2:nj) + g%vx(:, 1:nj-1)) + &
            g%lx_i * 0.5 * (g%p(:, 2:nj) + g%p(:, 1:nj-1))
    flux_j = mass_j * 0.5 * (g%vx(1:ni-1, :) + g%vx(2:ni, :)) + &
            g%lx_j * 0.5 * (g%p(1:ni-1, :) + g%p(2:ni, :))
    
end subroutine calculate_x_momentum_flux_second_order

subroutine calculate_y_momentum_flux_second_order(g, mass_i, mass_j, flux_i, flux_j)
    type(t_grid), intent(inout) :: g
    real, intent(in) :: mass_i(:,:), mass_j(:,:)
    real, intent(inout) :: flux_i(:,:), flux_j(:,:)
    integer :: ni, nj

    ni = g%ni
    nj = g%nj

    flux_i = mass_i * 0.5 * (g%vy(:, 2:nj) + g%vy(:, 1:nj-1)) + &
            g%ly_i * 0.5 * (g%p(:, 2:nj) + g%p(:, 1:nj-1))
    flux_j = mass_j * 0.5 * (g%vy(1:ni-1, :) + g%vy(2:ni, :)) + &
            g%ly_j * 0.5 * (g%p(1:ni-1, :) + g%p(2:ni, :))
    
end subroutine calculate_y_momentum_flux_second_order

end module flux_calculation