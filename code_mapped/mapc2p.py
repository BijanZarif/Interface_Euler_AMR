# RING INCLUSION MAPPING FROM COMPUTATIONAL TO PHYSICAL DOMAIN FUNCTION
# Mauricio del Razo's extension of the paper:
# D. A. Calhoun, C. Helzel, and R. J. LeVeque, Logically rectangular grids and finite volume
# methods for pdes in circular and spherical domains, SIAM review, (2008), pp. 723-752.
# Changes here have to also be included in mapc2p.f90

# Define some global parameters (also need to be adjusted in setrun.py)
rsqr0 = 0.039999  # Radius of rectangle in computational grid
rout0 = 0.015     # Inner radius of ring inclusion
rinn0 = 0.010     # Outer radius of ring inclusion

# Define mapped grid function
# NOTE it requires grids multiples of (x,y)=(20,10) to match interfaces exactly
def mapc2p(xc,yc):
    from math import sqrt
    from pylab import sign
    
    # save global parameters in local variables
    rsqr = rsqr0
    rinn = rinn0
    rout = rout0
    
    # Initialize physical domain vectors
    xp = 1.0*xc
    yp = 1.0*yc
    for i in range(len(xc)):
      for j in range(len(xc[0])):

        xci = xc[i][j]
        yci = yc[i][j]
        
        # Unaffected part of computational grid
        if (abs(xci) >= rsqr or abs(yci) >= rsqr):
            xp[i][j] = 1.0*xci
            yp[i][j] = 1.0*yci
        else: 
            #! Scale computational grid: map rsqr to unit square
            xci = xc[i][j]/rsqr
            yci = yc[i][j]/rsqr
            
            # Calculate D,R, center for inner, middle and outer region
            d = max(abs(xci),abs(yci),1e-10)
            d = min(0.999999,d)
            D = rsqr*d/sqrt(2)
            rat = rinn/rsqr
            rat2 = (rout/rsqr)
            # For inner most gird lines inside rinn
            if d <= rat: # or use 0.5 to make sure it aligns with 
                R = rinn
            # For middle grid lines between rinn and rout
            elif d <= rat2: # or use rat2
                R =  rsqr*d
            # For outer grid lines outside rout (note D has to be redefined)
            else:
                R = rout*((1.-rat2)/(1.-d))**(1.0/rat2 + 1.0/2.0) 
                D = rout/sqrt(2) + (d-rat2)*(rsqr-rout/sqrt(2))/(1.-rat2)
            center = D - sqrt(R**2 - D**2)
            
            # Do mapping (see paper) Have to do for the 3 sections from the diagonals
            xp[i][j] = (D/d)*abs(xci)
            yp[i][j] = (D/d)*abs(yci)
            if abs(xci) >= abs(yci):
                xp[i][j] =  center + sqrt(R**2 - yp[i][j]**2)
            if abs(yci) >= abs(xci):
                yp[i][j] = center + sqrt(R**2 - xp[i][j]**2)
            # Give corresponding sign
            xp[i][j] = sign(xci)*xp[i][j]
            yp[i][j] = sign(yci)*yp[i][j]
    
    # Uncomment to return no mapped grid at all        
    #xp = 1.0*xc
    #yp = 1.0*yc
    return xp,yp

# Function to test the mapping, mx and my are the number of grid cells in x and y
# Try calling test(40,20) from ipython after running run mapc2p.py.
def test(mx,my):
    import numpy
    from matplotlib import pyplot as plt
    from matplotlib import cm
    ## Test for unit
    sc = 100.0 # Scaling
    routS = rout0*sc
    rinnS = rinn0*sc
    xc_1d = numpy.linspace(-0.05, 0.05, mx+1)
    yc_1d = numpy.linspace(0, 0.05, my+1)
    XC,YC = numpy.meshgrid(xc_1d, yc_1d)
    XP,YP = mapc2p(XC,YC)
    XC = XC*sc
    YC = YC*sc
    XP = XP*sc
    YP = YP*sc
    eps = 1e-12
    half_dx = 0.5*sc*(xc_1d[1] - xc_1d[0])
    half_dy = 0.5*sc*(yc_1d[1] - yc_1d[0])
    yc2 = routS - 2*half_dy
    yc3 = routS
    xc3 = routS - 2*half_dx
    xc4 = routS
    Z = numpy.zeros(XC.shape)
    Z = numpy.where((YC>=yc2-eps) & (YC<=yc3-eps), 0.6, Z)
    Z2 = numpy.where((XC>xc3-eps) & (XC<xc4), 1., Z)

    plt.figure(1,figsize=(10,5))
    plt.clf()
    plt.plot(XP,YP,color=[.8,.8,.8])
    plt.plot(XP.T,YP.T,color=[.8,.8,.8])
    plt.pcolor(XP,YP,Z,cmap=cm.Greens)
    plt.pcolor(XP,YP,Z2,cmap=cm.Greens)
    plt.clim(0,2)
    
    # Plot interface
    xint = numpy.linspace(0,routS,10000)
    yint = numpy.sqrt(routS**2 - xint**2)
    xint2 = numpy.linspace(0,rinnS,10000)
    yint2 = numpy.sqrt(rinnS**2 - xint2**2)
    plt.plot(-xint,yint, '-k', linewidth=3, label='Interface 1')
    plt.plot(xint,yint, '-k', linewidth=3)
    plt.plot(-xint2,yint2, '--k', linewidth=2, label='Interface 2')
    plt.plot(xint2,yint2, '--k', linewidth=2)
    plt.title("Physical domain", fontsize=20)
    plt.axes().set_aspect('equal')
    plt.xlim([-0.05*sc,0.05*sc])
    plt.ylim([0.0,0.05*sc])
    plt.xlabel(r"$ x_p (cm)$", fontsize=18)
    plt.ylabel(r"$ y_p (cm)$", fontsize=18)
    plt.legend()
 
    # Plot figures
    plt.figure(2,figsize=(10,5))
    plt.clf()
    plt.plot(XC,YC,color=[.8,.8,.8])
    plt.plot(XC.T,YC.T,color=[.8,.8,.8])
    plt.pcolor(XC,YC,Z,cmap=cm.Greens)
    plt.pcolor(XC,YC,Z2,cmap=cm.Greens)
    plt.clim(0,2)
    xcart = [-routS, -routS, routS, routS]
    ycart = [0., routS, routS,  0.]
    xcart2 = [-rinnS, -rinnS, rinnS, rinnS]
    ycart2 = [0., rinnS, rinnS,  0.]
    plt.plot(xcart,ycart, '-k', linewidth=3, label='Interface 1')    
    plt.plot(xcart2,ycart2, '--k', linewidth=2, label='Interface 2')
    plt.title("Computational domain" , fontsize=20)
    plt.axes().set_aspect('equal')
    plt.xlim([-0.05*sc,0.05*sc])
    plt.ylim([0.0,0.05*sc])
    plt.xlabel(r"$ x_c$", fontsize=18)
    plt.ylabel(r"$ y_c$", fontsize=18)
    plt.legend()

    plt.show()

