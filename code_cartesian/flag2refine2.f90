! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Default version computes spatial difference dq in each direction and
! for each component of q and flags any point where this is greater than
! the tolerance tolsp.  This is consistent with what the routine errsp did in
! earlier versions of amrclaw (4.2 and before).
!
! This routine can be copied to an application directory and modified to
! implement some other desired refinement criterion.
!
! Points may also be flagged for refining based on a Richardson estimate
! of the error, obtained by comparing solutions on the current grid and a
! coarsened grid.  Points are flagged if the estimated error is larger than
! the parameter tol in amr2ez.data, provided flag_richardson is .true.,
! otherwise the coarsening and Richardson estimation is not performed!  
! Points are flagged via Richardson in a separate routine.
!
! Once points are flagged via this routine and/or Richardson, the subroutine
! flagregions is applied to check each point against the min_level and
! max_level of refinement specified in any "region" set by the user.
! So flags set here might be over-ruled by region constraints.
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)    
!
! tolsp = tolerance specified by user in input file amr2ez.data, used in default
!         version of this routine as a tolerance for spatial differences.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                            tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)

    use regions_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,j,m
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: dqi(meqn), dqj(meqn), dq(meqn)
    
    ! Varaibles required to calculate pressure @maojrs
    real(kind=8) :: uip, vip, ekip, gammaip, pinfip, pip
    real(kind=8) :: uim, vim, ekim, gammaim, pinfim, pim
    real(kind=8) :: dpi, dpj, dq2

    ! Initialize flags
    amrflags = DONTFLAG
    
    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    ! This information is not needed for the default flagging based on
    ! undivided differences, but might be needed in a user's version.
    ! Note that if you want to refine only in certain space-time regions,
    ! it may be easiest to use the "regions" feature.  The flags set here or
    ! in the Richardson error estimator are potentially modified by the
    ! min_level and max_level specified in any regions.

    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy
        
        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

            ! -----------------------------------------------------------------
                dq = 0.d0
                dqi = abs(q(:,i+1,j) - q(:,i-1,j))
                dqj = abs(q(:,i,j+1) - q(:,i,j-1))
                dq = max(dq,dqi,dqj)
                
                ! Calculate dpi (pressure difference in x) @maojrs
                uip = q(2,i+1,j)/q(1,i+1,j)
                vip = q(3,i+1,j)/q(1,i+1,j)
                ekip = 0.5*q(1,i+1,j)*(uip**2 + vip**2)
                gammaip = aux(1,i+1,j)
                pinfip = aux(2,i+1,j)
                pip = (q(4,i+1,j) - ekip)*(gammaip - 1.0) - gammaip*pinfip
                
                uim = q(2,i-1,j)/q(1,i-1,j)
                vim = q(3,i-1,j)/q(1,i-1,j)
                ekim = 0.5*q(1,i-1,j)*(uim**2 + vim**2)
                gammaim = aux(1,i-1,j)
                pinfim = aux(2,i-1,j)
                pim = (q(4,i-1,j) - ekim)*(gammaim - 1.0) - gammaim*pinfim
                
                dpi = abs(pip - pim)
                
                ! Calculate dpj (pressure difference in y) @maojrs
                uip = q(2,i,j+1)/q(1,i,j+1)
                vip = q(3,i,j+1)/q(1,i,j+1)
                ekip = 0.5*q(1,i,j+1)*(uip**2 + vip**2)
                gammaip = aux(1,i,j+1)
                pinfip = aux(2,i,j+1)
                pip = (q(4,i,j+1) - ekip)*(gammaip - 1.0) - gammaip*pinfip
                
                uim = q(2,i,j-1)/q(1,i,j-1)
                vim = q(3,i,j-1)/q(1,i,j-1)
                ekim = 0.5*q(1,i,j-1)*(uim**2 + vim**2)
                gammaim = aux(1,i,j-1)
                pinfim = aux(2,i,j-1)
                pim = (q(4,i,j-1) - ekim)*(gammaim - 1.0) - gammaim*pinfim
                
                dpj = abs(pip - pim)
                
                dq2 = 0.d0
                dq2 = max(dq2,dpi,dpj)
! 
!                 ! default checks all components of undivided difference:
!                 do m=1,meqn
!                     if (dq(m) > tolsp) then
!                         amrflags(i,j) = DOFLAG
!                         cycle x_loop
!                     endif
!                 enddo
                 
                ! Override default flagging by pressure pressure tolerance check @maojrs
                if (dq2 > tolsp) then
                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif

        enddo x_loop
    enddo y_loop
   
end subroutine flag2refine2
