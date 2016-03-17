c
c =========================================================
      subroutine rpn2(ixy,maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 2D Euler equations (normal solver) 
c     # using HLLC - Tamman-exact hybrid approximate Riemann solver
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
        
      implicit double precision (a-h,o-z)
      dimension   ql(meqn, 1-mbc:maxmx+mbc)
      dimension   qr(meqn, 1-mbc:maxmx+mbc)
      dimension   qm(meqn, 1-mbc:maxmx+mbc)
      dimension   qml(meqn, 1-mbc:maxmx+mbc)
      dimension   qmr(meqn, 1-mbc:maxmx+mbc)
      dimension   fluxl(meqn,1-mbc:maxmx+mbc)
      dimension   fluxr(meqn,1-mbc:maxmx+mbc)
      dimension    s(mwaves, 1-mbc:maxmx+mbc)
      dimension wave(meqn, mwaves, 1-mbc:maxmx+mbc)
      dimension amdq(meqn, 1-mbc:maxmx+mbc)
      dimension apdq(meqn, 1-mbc:maxmx+mbc)
      dimension auxl(maux,1-mbc:maxmx+mbc)
      dimension auxr(maux,1-mbc:maxmx+mbc)
      
      dimension q2l_rot(1-mbc:maxmx+mbc),q3l_rot(1-mbc:maxmx+mbc)
      dimension q2r_rot(1-mbc:maxmx+mbc),q3r_rot(1-mbc:maxmx+mbc)
      dimension grid_nx(1-mbc:maxmx+mbc),grid_ny(1-mbc:maxmx+mbc)
      

c     # local storage
c     ---------------
!       parameter (max2 = 20002)  !# assumes at most 2000 grid points with mbc=2
     
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
      common /param2/ dt2, dx2,dy2
c
c     # Dimensional splitting
      if (ixy.eq.1) then
          mu = 2
          mv = 3
          map = 8
          isig = 1
          dxx = dx2
      else
          mu = 3
          mv = 2
          isig = -1
          map = 11
          dxx = dy2
      endif
      
      do i=2-mbc,mx+mbc
      	  grid_nx(i) = auxl(map,i)
	  grid_ny(i) = auxl(map+1,i)
	  q2l_rot(i) = ql(2,i)*grid_nx(i) + ql(3,i)*grid_ny(i)
	  q2r_rot(i-1) = qr(2,i-1)*grid_nx(i) + qr(3,i-1)*grid_ny(i)
	  q3l_rot(i) = -ql(2,i)*grid_ny(i) + ql(3,i)*grid_nx(i)
	  q3r_rot(i-1) = -qr(2,i-1)*grid_ny(i) + qr(3,i-1)*grid_nx(i)
      end do
      
c     # Iterate over line of grid cells to solve Riemann problem
      do 20 i=2-mbc,mx+mbc
         ! Calculate current cell
         xcell = -10.0 + (i-0.5d0)*.052
         if (ixy .eq. 2) then
            xcell = (i-0.5d0)*.05
         end if
         
         !Obtain SGEOS parameters from aux arrays
         gammal = auxr(1,i-1)
         gammar = auxl(1,i)
         gamma1l = gammal - 1.0
         gamma1r = gammar - 1.0
         pinfl = auxr(2,i-1)
         pinfr = auxl(2,i)
         omel = auxr(3,i-1)
         omer = auxl(3,i)
         
         ! Compute main quantities
         ! Densities
         rho_l = qr(1,i-1)
         rho_r = ql(1,i)
	  
         ! Velocities 
         ul = q2r_rot(i-1)/rho_l
         ur = q2l_rot(i)/rho_r
         vl = q3r_rot(i-1)/rho_l
         vr = q3l_rot(i)/rho_r
         
         ! Kinetic Energy
         ek_l = 0.5*rho_l*(ul**2 + vl**2)
         ek_r = 0.5*rho_r*(ur**2 + vr**2)
         ! Pressures (Use Tait EOS on water and/or plastic, SGEOS on air or SGEOS on both)
         pl = gamma1l*(qr(4,i-1) - ek_l) 
         pl = pl/(1.0 - omel*rho_l) - pinfl*gammal
         pr = gamma1r*(ql(4,i) - ek_r) 
         pr = pr/(1.0 - omer*rho_r) - pinfr*gammar
    
         ! Compute left and right speeds
         cl = dsqrt(gammal*(pl + pinfl)/rho_l)
         cr = dsqrt(gammar*(pr + pinfr)/rho_r)

         ! Compute the speed of left and right HLLC wave
         Sl = min(ul - cl,ur - cr) ! u(i) - a(i)
         Sr = max(ul + cl,ur + cr) ! u(i) + a(i),        
         s(1,i) = 1.d0*Sl
         s(3,i) = 1.d0*Sr
         
         ! Compute HLLC middle speed state (see research notebook)
         Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
         Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
         s(2,i) = 1.d0*Sm
         
         ! Calculate solution values in middle states
         do j=1,meqn
             qml(j,i) = rho_l*(Sl - ul)/(Sl - Sm)
             qmr(j,i) = rho_r*(Sr - ur)/(Sr - Sm)
         end do
         ! Add multiplicative terms to momentum ones
         qml(2,i) = Sm*qml(2,i)
         qmr(2,i) = Sm*qmr(2,i)
         qml(3,i) = vl*qml(3,i)
         qmr(3,i) = vr*qmr(3,i)
         ! Add second terms to energy one (see Toro pg. 325)
         qml(4,i) = qml(4,i)*(qr(4,i-1)/rho_l + 
     & (Sm - ul)*(Sm + pl/(rho_l*(Sl - ul))))
         qmr(4,i) = qmr(4,i)*(ql(4,i)/rho_r + 
     & (Sm - ur)*(Sm + pr/(rho_r*(Sr - ur))))
	! Ends HLLC solver
         
         ! Run alternative Riemann process if we are at interface
         if ((gammal*gammar .eq. gammagas*gammawat)) then !
            ! Do exact Tamman Riemann solver at interface and in water
            ! Newton's ,method to obtain correct pressure and velocity pstar,ustar
            pstar = 0.5*(pl + pr)
            pold = pstar + 10
            do while (abs(pstar - pold) > 0.0001)
              pold = pstar
              ! Call function to solve exact riemann solver for Tamman EOS
              CALL phi_exact(gammal,gammar,pr,pl,rho_r,rho_l,ul,ur,
     & pinfl,pinfr,pstar,phi,phi_prime,rhos_l,rhos_r,ustar)
              pstar = pstar - phi/phi_prime
            end do
            
            ! Compute the speed of left right and middle wave (See MJ IVINGS paper)
            betal = (pl + pinfl)*(gammal - 1.0)/(gammal + 1.0)
            betar = (pr + pinfr)*(gammar - 1.0)/(gammar + 1.0)
            alphal = 2.0/(rho_l*(gammal + 1.0))
            alphar = 2.0/(rho_r*(gammar + 1.0))
            Sm = 1.0*ustar
            Sl = ul - dsqrt((pstar + pinfl + betal)/alphal)/rho_l
            Sr = ur + dsqrt((pstar + pinfr + betar)/alphar)/rho_r
            
            ! Always move to a Lagrangian frame of reference at interface
            s(1,i) = 1.d0*Sl - 1.0*ustar
            s(2,i) = 0.0*ustar
            s(3,i) = 1.d0*Sr - 1.0*ustar
            
            ! Calculate middle states using pstar and ustar of exact solution
            bl = (gammal + 1.0)/(gammal - 1.0)
            br = (gammar + 1.0)/(gammar - 1.0)
            ! Calculate densities, momentums and energys
            qml(1,i) = rhos_l !rho_l*(1 + bl*pstar/pl)/(pstar/pl + bl)
            qmr(1,i) = rhos_r !rho_r*(1 + br*pstar/pr)/(pstar/pr + br)
            qml(2,i) = qml(1,i)*ustar
            qmr(2,i) = qmr(1,i)*ustar
            qml(3,i) = qml(1,i)*vl !*rho_l*(Sl - ul)/(Sl - ustar)
            qmr(3,i) = qmr(1,i)*vr!*rho_r*(Sr - ur)/(Sr - ustar)
            qml(4,i) = (pstar + gammal*pinfl)/(gammal - 1.0) + 
     & 0.5*(qml(2,i)**2 + qml(3,i)**2)/qml(1,i)
            qmr(4,i) = (pstar + gammar*pinfr)/(gammar - 1.0) + 
     & 0.5*(qmr(2,i)**2 + qmr(3,i)**2)/qmr(1,i)
         end if
         ! Ends exact solver
 
c        # Compute the 3 waves with values obtained from HLLC or exact solver.
c        j index over q variables
         do j=1,meqn
             q_l = qr(j,i-1)
             q_r = ql(j,i)
             if (j .eq. 2) then
	       q_l = q2r_rot(i-1)
               q_r = q2l_rot(i)
             else if (j .eq. 3) then
	       q_l = q3r_rot(i-1)
               q_r = q3l_rot(i)
             end if
             wave(j,1,i) = qml(j,i) - q_l
             wave(j,2,i) = qmr(j,i) - qml(j,i)
             wave(j,3,i) = q_r - qmr(j,i) 
         end do
	
         
   20    continue
   
c   ! Rotate waves back to physical domain to use the mapped grid
	do i=2-mbc, mx+mbc
	    do mw=1,mwaves
	      wave2 = wave(2,mw,i)*grid_nx(i) - wave(3,mw,i)*grid_ny(i)
	      wave3 = wave(2,mw,i)*grid_ny(i) + wave(3,mw,i)*grid_nx(i)

	      wave(2,mw,i) = 1.0*wave2
	      wave(3,mw,i) = 1.0*wave3
	      
	      ! Scales speed by relative length of edge of mapped grid	
	      s(mw,i) = s(mw,i)*auxr(map+2,i)
	      end do
	  end do
      
c     ! Compute the leftgoing and rightgoing flux differences
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
               if (s(mw,i) .lt. 0.d0) then
                   amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
   90          continue
  100       continue
      go to 900
c
c-----------------------------------------------------
  900 continue
      return
      end

