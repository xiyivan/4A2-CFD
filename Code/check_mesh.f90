    
    subroutine check_mesh(g)

!   Check the cell area and facet length calculations before attempting to
!   solve the flow, make sure you do this for both the "bump" and "bend" test
!   cases

!   Explicitly declare the required variables
    use types
    implicit none
    type(t_grid), intent(inout) :: g
    real :: area_min, dx_error, dy_error, tol
    integer :: ni, nj

!   Get the size of the mesh and store locally for convenience
    ni = g%ni; nj = g%nj;

!   Exact checking of floating point numbers never goes well, define a
!   small tolerance to use for a comparative operation instead
    tol = 1e-4 * g%l_min

!   Check that all of the cell areas are positive, either with the intrinsic
!   "minval" function or with nested do loops. Print the output to the screen
!   and flag negative numbers as an error with an if statement to "stop" the
!   program
!   INSERT

    area_min = minval(g%area)
    

    write(6,*) "Minimum area ", area_min

    if (area_min < 0) then
        stop "None positive area detected, program exiting"
    else 
        write(6,*) "Area check passed"
    end if

!   END INSERT 
        

!   Next check that the sum of the edge vectors around every quadrilateral is 
!   very nearly zero in both the x and y-coordinate directions. You can
!   complete this with some elementwise addition of the arrays and use of the
!   "maxval" and "abs" intrinsic functions.
!   INSERT

    dx_error = maxval(abs(g%lx_i(1:ni-1, :) - g%lx_i(2:ni, :) + &
    g%lx_j(:, 1:nj-1) - g%lx_j(:, 2:nj)))

    dy_error = maxval(abs(g%ly_i(1:ni-1, :) - g%ly_i(2:ni, :) + &
    g%ly_j(:, 1:nj-1) - g%ly_j(:, 2:nj)))

    write(6,*) "minimum round-loop vecor error as below"
    write(6,*) "dx_error: ", dx_error
    write(6,*) "dy_error: ", dy_error

    if (dx_error > tol .OR. dy_error > tol) then
        stop "Large round-loop vector error detected, program exiting"
    else 
        write(6,*) "Loop check passed"
    end if

!   END INSERT

!   It may be worthwhile to complete some other checks, the prevous call to
!   the "write_output" subroutine has written a file that you can read and
!   postprocess using the Python script plot_mesh.py. This program also has
!   access to all of the mesh parameters used within the solver that you could
!   inspect graphically.

!   Print a blank line
    write(6,*)

    end subroutine check_mesh
