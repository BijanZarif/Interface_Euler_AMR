c     ============================================
      subroutine b4step2(mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called from claw2 before each call to step2.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.
 
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      
      common /param2/ dt2, dx2
      
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ rhog,rhop,rhow
      
      dt2 = 1.d0*dt
      dx2 = 1.d0*dx
c
	  ! Write q variable to aux 4,5,6 and 7 array
      ! This will be required by the riemann transverse 
      ! solver (rpt2acoustic in rp folder)
      do ii=1-mbc,mx+mbc
		do jj=1-mbc,my+mbc
		  aux(4,ii,jj) = q(1,ii,jj)
		  aux(5,ii,jj) = q(2,ii,jj)
		  aux(6,ii,jj) = q(3,ii,jj)
		  aux(7,ii,jj) = q(4,ii,jj)
		end do
      end do 
      
c
      return
      end

