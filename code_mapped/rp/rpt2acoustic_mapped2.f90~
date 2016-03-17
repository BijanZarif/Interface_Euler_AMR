! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================

      implicit double precision (a-h,o-z)
!
!     # Riemann solver in the transverse direction for the Euler equations
!     # with Tammann EOS using acoustic approximation to deal with interfaces.
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
!
!     # Uses approximation to acoustic eqs to calculate downward and upward
!     # contributions with impedance and sound speed where there is an interface
!
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension aux1(maux, 1-mbc:maxm+mbc)
      dimension aux2(maux, 1-mbc:maxm+mbc)
      dimension aux3(maux, 1-mbc:maxm+mbc)
!
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
      dimension waveb(4,3),sb(3)
      parameter (maxm2 = 1800)  !# assumes at most 200x200 grid with mbc=2
!
      double precision u2v2(3,-6:maxm2+7), &
            u(3,-6:maxm2+7), v(3,-6:maxm2+7), &
            enth(3,-6:maxm2+7),a(3,-6:maxm2+7), &
            g1a2(3,-6:maxm2+7),euv(3,-6:maxm2+7)

      if (mbc > 7 .OR. maxm2 < maxm) then
         write(6,*) 'need to increase maxm2 or 7 in rpn'
         stop
      endif
!
      gamma1 = gammagas - 1.d0

      if (ixy.eq.1) then
          mu = 2
          mv = 3
          map = 11
        else
          mu = 3
          mv = 2
          map = 8
        endif
        

!     Calculate relevant quantities in the three transverse cells 3-up, 2-middle, 1 down 
      do 10 i = 2-mbc, mx+mbc

         if (imp.eq.1) then 
!            # asdq = amdq, moving to left
             i1 = i-1
         else
!            # asdq = apdq, moving to right
             i1 = i
         endif

        ! Obtain current data from aux arrays aux1,aux2,aux3
        ! 3-up 2-middle 1-down
        
	  ! Densities
        rho3 = aux3(4,i1)
        rho2 = aux2(4,i1)
        rho1 = aux1(4,i1)
        
	  ! Normal velocities
        u(3,i1) = aux3(5,i1)/rho3
        u(2,i1) = aux2(5,i1)/rho2
        u(1,i1) = aux1(5,i1)/rho1
	  
	  ! Transverse velocities
        v(3,i1) = aux3(6,i1)/rho3
        v(2,i1) = aux2(6,i1)/rho2
        v(1,i1) = aux1(6,i1)/rho1
        
      ! Internal energies
	    ene3 = aux3(7,i1) 
        ene2 = aux2(7,i1) 
        ene1 = aux1(7,i1) 
        
	  ! Kinetc energies 
	    enek3 = 0.5*(u(3,i1)**2 + v(3,i1)**2)/rho3
        enek2 = 0.5*(u(2,i1)**2 + v(2,i1)**2)/rho2
        enek1 = 0.5*(u(1,i1)**2 + v(1,i1)**2)/rho1

	 ! Obtain parameters from aux array
		gamma3 = aux3(1,i1)
		gamma2 = aux2(1,i1)
		gamma1 = aux1(1,i1)
		pinf3 = aux3(2,i1)
		pinf2 = aux2(2,i1)
		pinf1 = aux1(2,i1)
	 
	 ! Calculate pressures
	p3 = (gamma3 - 1)*(ene3 - enek3) - gamma3*pinf3
	p2 = (gamma2 - 1)*(ene2 - enek2) - gamma2*pinf2
	p1 = (gamma1 - 1)*(ene1 - enek1) - gamma1*pinf1
	
	 ! Calculate enthalpies
	enth(3,i1) = (ene3 + p3)/rho3 
	enth(2,i1) = (ene2 + p2)/rho2
	enth(1,i1) = (ene1 + p1)/rho1
	
	 ! Calculate u^2 + v^2 array
	u2v2(3,i1) = u(3,i1)**2 + v(3,i1)**2
	u2v2(2,i1) = u(2,i1)**2 + v(2,i1)**2
	u2v2(1,i1) = u(1,i1)**2 + v(1,i1)**2
	
	! Calculate speeds c
	a(3,i1) = dsqrt(gamma3*(p3 + pinf3)/rho3)
	a(2,i1) = dsqrt(gamma2*(p2 + pinf2)/rho2)
	a(1,i1) = dsqrt(gamma1*(p1 + pinf1)/rho1)
	
	 ! Calculate (gamma-1)/c^2 array
	g1a2(3,i1) = (gamma3 -1.0)/a(3,i1)**2	 
	g1a2(2,i1) = (gamma2 -1.0)/a(2,i1)**2
	g1a2(1,i1) = (gamma1 -1.0)/a(1,i1)**2

	 ! Calculate H-u^2-v^2 array 
	euv(3,i1) = enth(3,i1) - u2v2(3,i1) 
	euv(2,i1) = enth(2,i1) - u2v2(2,i1)
	euv(1,i1) = enth(1,i1) - u2v2(1,i1)

	
   10    continue
!
         ! Do transverse solver for upward direction using Bi,j+1 matrix
         ! And for downward direction Bi,j-1 matrix, no Roe averaging done
         do 20 i = 2-mbc, mx+mbc

		 if (imp.eq.1) then 
			 ! # asdq = amdq, moving to left
             i1 = i-1
         else
			 ! # asdq = apdq, moving to right
             i1 = i
         endif

		 !Define direction of normal to grid edge normals for downgoing fluctuation
	 	nxm = aux2(map,i1)
	 	nym = aux2(map+1,i1) 
	 
		!Define direction of normal to grid edge normals for upgoing fluctuation
	 	nxp = aux3(map,i1)
	 	nyp = aux3(map+1,i1) 

!        # jumps in asdq:
         drho   = asdq(1,i1)
         du     = asdq(2,i1)
         dv     = asdq(3,i1)
	    
	    cond = 1
	    ! Check if we are on an Horizontal interface in transverse direction (and inside water)
	    if ((gamma1*gamma2 .eq. gammagas*gammawat) .or. (gamma2*gamma3 .eq. gammagas*gammawat) & 
	       .or. (gamma1*gamma2 .eq. gammawat*gammawat) .or. (gamma2*gamma3 .eq. gammawat*gammawat)) then
	      ! Do normal transverse ROE solver
	      cond = 0 ! 1 = uses traditional solver for air
	    end if
	    
	    ! Check if we are at corner of interface to not do transverse solver (avoids instabilities)
	    if (aux1(15,i1) + aux2(15,i1) + aux3(15,i1) > 0.5) then
            cond = 2
        end if
	    
	    ! Do transverse solvers for interface case (acoustic approximation as in paper)
	    ! Note velocity of "contact discontinuity" substracted (- v({3,1},i))
        if (cond .eq. 0) then 
			cd = a(1,i1)
			cu = a(3,i1)			
	            
	  !        # transmitted part of down-going wave:
	      a1 = (cd*drho + nxm*du + nym*dv) / &
	      (cd + cu)

	  !        # transmitted part of up-going wave:
	      a23 = (cu*drho -nxp*du - nyp*dv) / &
	      (cd + cu)
  !
  !           # Use acoustics eqs to calculate acoustic waves;
  !           # Contact discontinuity and shear wave neglected 
  !           # (Consider we are using Lagrangian frame of ref, since grid is fixed)
  !           # compute the flux differences bmasdq (Wave x /pm speed a)
	      cdscale = cd * aux2(map+2,i1) ! Scales speed by relative length of edge of mapped grid
          bmasdq(1,i1) = cdscale * a1 ! Acoustic version of aocustic wave
	      bmasdq(2,i1) = -cdscale * cd * nxm * a1
	      bmasdq(3,i1) = -cdscale * cd * nym * a1  ! Acoustic version of aocustic wave
	      bmasdq(4,i1) = 0.d0
	      
  !           # compute the flux differences bpasdq
          cuscale = cu*aux3(map+2,i1)
	      bpasdq(1,i1) = cuscale * a23 ! Acoustic version of aocustic wave
	      bpasdq(2,i1) = cuscale * cu * a23 * nxp
	      bpasdq(3,i1) = cuscale * cu * a23 * nyp ! Acoustic version of aocustic wave
	      bpasdq(4,i1) = 0.d0
    
	    end if
	    
	    ! Do usual no-interface transverse solver
	    !(without ROE avgs (there could be an interface in normal direction)
	    on = 1
	    if ((cond .eq. 1) .or. (cond .eq. 2)) then
          if (cond .eq. 2) then
            on = 0
          end if
	      a3 = g1a2(2,i) * (euv(2,i)*asdq(1,i) &
		    + u(2,i)*asdq(mu,i) + v(2,i)*asdq(mv,i) - asdq(4,i))
	      a2 = asdq(mu,i) - u(2,i)*asdq(1,i)
	      a4 = (asdq(mv,i) + (a(2,i)-v(2,i))*asdq(1,i) - a(2,i)*a3) &
		    / (2.d0*a(2,i))
	      a1 = asdq(1,i) - a3 - a4
  !
	      waveb(1,1) = on*a1
	      waveb(mu,1) = on*a1*u(2,i)
	      waveb(mv,1) = on*a1*(v(2,i)-a(2,i))
	      waveb(4,1) = on*a1*(enth(2,i) - v(2,i)*a(2,i))
	      sb(1) = on*(v(2,i) - a(2,i))
  !
	      waveb(1,2) = on*a3
	      waveb(mu,2) = on*(a3*u(2,i) + a2)
	      waveb(mv,2) = on*a3*v(2,i)
	      waveb(4,2) = on*(a3*0.5d0*u2v2(2,i) + a2*u(2,i))
	      sb(2) = on*v(2,i)
  !
	      waveb(1,3) = on*a4
	      waveb(mu,3) = on*a4*u(2,i)
	      waveb(mv,3) = on*a4*(v(2,i)+a(2,i))
	      waveb(4,3) = on*a4*(enth(2,i)+v(2,i)*a(2,i))
	      sb(3) = on*(v(2,i) + a(2,i))
      end if
  !      
  !           # compute the flux differences bmasdq and bpasdq
  !
	  
	  do 50 m=1,meqn
		bmasdq(m,i) = 0.d0
		bpasdq(m,i) = 0.d0
		do 50 mw=1,3
		    bmasdq(m,i) = bmasdq(m,i) &
				+ dmin1(sb(mw), 0.d0) * waveb(m,mw)
		    bpasdq(m,i) = bpasdq(m,i) &
				+ dmax1(sb(mw), 0.d0) * waveb(m,mw)
    50             continue
!                 
   20       continue
!
      return
      end
