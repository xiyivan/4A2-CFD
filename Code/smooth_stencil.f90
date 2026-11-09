
    module smooth_stencil

!   Packaging a subroutine in a module allows it to recieve the data
!   conveniently as assumed shape arrays
    
    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine smooth_array(av,prop,corr)

!   This subroutine smooths "prop" to stabilise the calculation, the basic 
!   solver uses second order smoothing, many improvements are possible.

!   Explicitly declare the required variables
    use types
    implicit none
    type(t_appvars), intent(in) :: av
    real, intent(inout) :: prop(:,:), corr(:,:)
    real, dimension(size(prop,1),size(prop,2)) :: prop_avg, corr_total
    integer :: ni, nj, i, j

!   Get the block size and store locally for convenience
    ni = size(prop,1); nj = size(prop,2)

!   Calculate the average values at the nodes in the interior region of the
!   mesh, use the four neighbouring nodes in the plus and minus i and 
!   j-directions.
!   INSERT

    do j = 2, nj-1
        do i = 2, ni-1
            prop_avg(i, j) = (prop(i+1, j) + &
                        prop(i-1, j) + &
                        prop(i, j+1) + &
                        prop(i, j-1))/4
        end do
    end do
    
!   END INSERT

!   Edge values are also averaged in both the i and j-directions. Parallel to
!   the boundary the averaging is centred, the averages of two nodes are taken
!   either side of the current point. Perpendicular to the boundary the
!   algorithm is one-sided, the value at the current point is extrapolated
!   from the values at two nodes away from the boundary point.
!   INSERT

    do i = 2, ni-1
        prop_avg(i, 1) = (prop(i+1, 1) + prop(i-1, 1) + 2*prop(i, 2) - prop(i,3))/3
        prop_avg(i, nj) = (prop(i+1, nj) + prop(i-1, nj) + 2*prop(i, nj-1) - prop(i,nj-2))/3
    end do

    do j = 2, nj-1
        prop_avg(1, j) = (prop(1, j+1) + prop(1, j-1) + 2*prop(2, j) - prop(3, j))/3
        prop_avg(ni, j) = (prop(ni, j+1) + prop(ni, j-1) + 2*prop(ni-1, j) - prop(ni-2, j))/3
    end do

!   END INSERT

!   The corner values are not currently smoothed
    prop_avg([1,ni],[1,nj]) = prop([1,ni],[1,nj])

!   use correction to reduce artificial viscosity
    if (av%fcorr == 0.0) then
        corr = 0.0
    else
        corr_total = av%fcorr * (prop-prop_avg)
        corr = 0.99 * corr + 0.01 * corr_total
    end if

!   Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!   take (1-sfac) * the calculated value of the property + sfac * the average 
!   of the surrounding values. 
!   INSERT

    prop = (1-av%sfac) * prop + av%sfac * (prop_avg + corr)

!   END INSERT

    end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module smooth_stencil


