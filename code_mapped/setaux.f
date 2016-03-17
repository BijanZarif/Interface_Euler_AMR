c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays to save Tammann EOS parameters, mapped grid normals,
c     # capacity function and location of corners.

      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      
      ! Arrays to temporarily store computational and physical corners of grid cells
      dimension xccorn(4),yccorn(4),xpcorn(4),ypcorn(4)
      double precision norm
      
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
      common /mappedparam/rsqr,rout,rinn

c     #   aux(1,i,j) = gamma in SGEOS ot Tamman EOS p = rho*e*(gamma - 1) - gamma*P_inf
c     #   aux(2,i,j) = P_inf in SGEOS ot Tamman EOS p = rho*e*(gamma - 1) - gamma*P_inf
c     #   aux(3,i,j) = Extra free parameter omega that can be used for improved EOS
                       !e.g. Polystyrene EOS based on Van der Waals approximation in Spender and Gilmore paper
      
c     #   These are defined in b4step.f
c     #   aux(4,i,j) saves q(1,i,j) into array for transverse solverse, (see b4step)
c     #   aux(5,i,j) saves q(2,i,j) into array for transverse solverse, (see b4step)
c     #   aux(6,i,j) saves q(3,i,j) into array for transverse solverse, (see b4step)
c     #   aux(7,i,j) saves q(4,i,j) into array for transverse solverse, (see b4step)

      
c     #   aux(8,i,j) is x component of normal at "left" boundary of grid point (i,j)
c     #   aux(9,i,j) is y component of normal at "left" boundary of grid point (i,j)
c     #   aux(10,i,j) =  is ratio of left edge to dy

c     #   aux(11,i,j) is x component of normal at "bottom" boundary of grid point (i,j)
c     #   aux(12,i,j) is y component of normal at "bottom" boundary of grid point (i,j)
c     #   aux(13,i,j) =  is ratio of bottom edge to dx

c     #   aux(14,i,j) = kappa  is ratio of cell area to dx*dx
c     #   aux(15,i,j) = 1 or zero. It is used to identify interface corners with a 1.
      
      
      ! Calculate square radius where the interface should be
      rad = rout 
      radin = rinn
      
      ! Loop over all cells
      do j=1-mbc,my + mbc 
		ycell = ylower + (j-0.5d0)*dy
		do i=1-mbc,mx + mbc 
	  		xcell = xlower + (i-0.5d0)*dx
            !============================================
	    	! CREATE INTERFACE USING DIFFERENT PARAMETERS FOR SGEOS (defined in setrun.py)
            if ((abs(xcell+0.0).le.rad) .and.(ycell.le.rad)) then
                aux(1,i,j) = gammawat
                aux(2,i,j) = pinfwat
                aux(3,i,j) = omewat
			! Skip plastic interface as in paper, it can be included if required.
   		    !else if ((xcell > radin).and.(ycell.le.radin)) then
            !      aux(1,i,j) = gammaplas
            !      aux(2,i,j) = pinfplas
            !      aux(3,i,j) = omeplas
            else 
                aux(1,i,j) = gammagas
                aux(2,i,j) = pinfgas
                aux(3,i,j) = omegas
            end if
            
            ! ============================================
            !BEGINS CALCULATING AND SAVING NORMALS, LENGTH RATIOS AND AREA RATIO FOR MAPPED GRID
            
            ! Computes corners of grid cell
c           # lower left corner:
	   	    xccorn(1) = xlower + float(i-1)*dx
			yccorn(1) = ylower + float(j-1)*dy
            call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1),
     &                  rsqr,rinn,rout)

c           # upper left corner:
		    xccorn(2) = xccorn(1)
		    yccorn(2) = yccorn(1) + dy
            call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2),
     &                  rsqr,rinn,rout)
c
c           # upper right corner:
		    xccorn(3) = xccorn(1) + dx
		    yccorn(3) = yccorn(1) + dy
            call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3),
     &                  rsqr,rinn,rout)
c
c           # lower right corner:
		    xccorn(4) = xccorn(1) + dx
		    yccorn(4) = yccorn(1)
            call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4),
     &                  rsqr,rinn,rout)
    	    
    	    ! Compute inner normals
    	    !left
    	    xn = ypcorn(2) - ypcorn(1)
    	    yn = -(xpcorn(2) - xpcorn(1))
    	    norm = dsqrt(xn**2 + yn**2)
    	    aux(8,i,j) = xn/norm
    	    aux(9,i,j) = yn/norm
    	    aux(10,i,j) = norm/dy
	    	!bottom
    	    xn = -(ypcorn(4) - ypcorn(1))
    	    yn = xpcorn(4) - xpcorn(1)
    	    norm = dsqrt(xn**2 + yn**2)
    	    aux(11,i,j) = xn/norm
    	    aux(12,i,j) = yn/norm
    	    aux(13,i,j) = norm/dx
    	    
    	    ! Computes area of grid cell using cross product of diagonal
    	    areap = (xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) 
    	    areap = areap - (ypcorn(3)-ypcorn(1))*(xpcorn(2)-xpcorn(4))
    	    areap = 0.5*abs(areap)
    	    aux(14,i,j) = areap/(dx*dy)
    	    
    	    
        end do
      end do
      
      ! Find corners cells and mark them as one all the others cells mark them as zero
      do j=2,my-1 
        do i=2,mx-1
            aux(15,i,j) = 0
            if ((aux(1,i,j)*aux(1,i+1,j)) .eq. (gammagas*gammawat)) then
              if ((aux(1,i+1,j+1)*aux(1,i+1,j+1)) .eq. 
     &                            (gammagas*gammagas)) then
                    aux(15,i,j) = 1
                    aux(15,i+1,j) = 1
                    aux(15,i,j+1) = 1
                    aux(15,i+1,j+1) = 1
              end if 
            end if
        end do
      end do
      
      return
      end


