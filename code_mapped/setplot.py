
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

# Note: To change plotted time scale, edit frametools.py in visclaw

import os
import numpy as np
from mapc2p import mapc2p
from matplotlib import rc
rc('text', usetex=True)


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Plot outline of interface withut mapping
    def aa(current_data):
      from pylab import linspace,plot,annotate,text
      from pylab import title, xlabel, ylabel, xticks, yticks, colorbar
      # Plot interface
      rout = 0.015
      rinn = 0.010
      x = linspace(-rout,rout, 100)
      y = np.sqrt(rout**2 - x**2)
      plot(x,y,'k',linewidth=3.0)
      #xinn = linspace(-rinn,rinn, 100)
      #yinn = np.sqrt(rinn**2 - xinn**2)
      #plot(xinn,yinn,'k',linewidth=3.0)
      # Chage title
      t = current_data.t
      tmicros = 1000000*t 
      title(r"AMR grids at time t = %10.2f $\mu s$" % tmicros, fontsize=16)
      # Change axes
      xlabel(r"$cm$", fontsize='16')
      ylabel(r"$cm$", fontsize='16')
      # Change ticks on axes (WATCHOUT IF DOMAIN OF SIMULATION IS CHANGED)
      xxticks = np.arange(-0.05, 0.05, 0.00999)
      labelsx = range(xxticks.size) 
      labelsx[:] = [x - 5 for x in labelsx]
      xticks(xxticks, labelsx)
      yyticks = np.arange(0.0, 0.03, 0.00999)
      labelsy = range(yyticks.size) 
      labelsy[:] = [y for y in labelsy]
      yticks(yyticks, labelsy)
      
    # Plot outline of interface
    def aa1DPSIcm(current_data):
      from pylab import linspace,plot,annotate,text,xlabel,ylabel
      #gcs = 2.0/200.0
      x = [-1.5,-1.5,1.5,1.5] 
      y = [-100,100,100,-100]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)
      xlabel('cm',fontsize='16')
      ylabel('psi',fontsize='16')
      xcav = [-3.0,3.0]
      ycav = [-14.334351113,-14.334351113] #Water vapour pressure for cavitation at room temp in 1atm=0 ref system
      plot(xcav,ycav,'b--')
      #plot(-8.0, 180000, 'vk', markersize=10) 
      #plot(-2.0, 180000, 'vk', markersize=10) 
      #plot(0.0, 180000, 'vk', markersize=10) 
      #plot(2.0, 180000, 'vk', markersize=10)
      text(-0.75,27,'Water',fontweight='bold',fontsize=20)
      #text(-0.8,285000,'PS',fontweight='bold',fontsize=20)
      text(-2.9,27,'Air',fontweight='bold',fontsize=20)
      text(1.6,27,'Air',fontweight='bold',fontsize=20)
      text(-1.45,-13,'Vapor pressure',fontsize=15,color='blue')
      
    # Function to calculate pressure when using Tammann EOS
    def Pressure(current_data):
        q = current_data.q   # solution when this function called
        aux = current_data.aux
        gamma = aux[0,:,:]
        gamma1 = aux[0,:,:] - 1.0
        pinf = aux[1,:,:]
        omega = aux[2,:,:]
        rho = q[0,:,:]           # density
        momx = q[1,:,:]          # momentum x
        momy = q[2,:,:]          # momentum y
        ene = q[3,:,:]           # energy
        P = gamma1*(ene - 0.5*(momx*momx + momy*momy)/rho)
        P = P - gamma*pinf
        return P

    # Figure for Density
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03] #'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Density'
    #plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    #plotitem.pcolor_cmin = 0.8
    #plotitem.pcolor_cmax = 3.0
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 1.0
    plotitem.pcolor_cmax = 2.0
    plotitem.show = True       # show on plot?
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    
    plotaxes.afteraxes = aa
    
    # Figure for momentum x
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Momentum x', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03] #'auto'
    plotaxes.ylimits = [-0.05,0.05] #'auto'
    plotaxes.title = 'Momentum x'
    #plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 160.0
    plotitem.show = True       # show on plot?
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    
    plotaxes.afteraxes = aa
    
    # Figure for momentum y
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Momentum y', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03]#'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Momentum y'
    #plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 160.0
    plotitem.show = True       # show on plot?
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    
    plotaxes.afteraxes = aa
    
    # Figure for Energy
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Energy', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.03,0.03]#'auto'
    plotaxes.ylimits = [-0.05,0.05]#'auto'
    plotaxes.title = 'Energy'
    #plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = 200000
    plotitem.pcolor_cmax = 400000
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    
    plotaxes.afteraxes = aa
    
    # Figure for Pressure
    # -------------------
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=4)
    plotfigure.kwargs = {'figsize':[8,3.7], 'tight_layout':True}
    #plotfigure.kwargs = {'figsize':[8,8], 'tight_layout':True} # For colorbar output

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('Pressure')
    plotaxes.xlimits = [-0.04,0.04] 
    plotaxes.ylimits = [0.001,0.035]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 90000
    plotitem.pcolor_cmax = 230000
    #plotitem.pcolor_cmap = colormaps.white_blue
    #white_green_cmap = colormaps.make_colormap({0.:'w', 0.35: '#54ED96', 0.7: '#31BCBC', 1.:'#005F8B'}) #5CDAE3
    #white_green_cmap = colormaps.make_colormap({0.:'w', 0.35: '#60E9D0', 0.7: '#3174B7', 1.:'#0B357F'}) #5CDAE3 
    white_green_cmap = colormaps.make_colormap({0.:'w', 0.35: '#AAFFEF', 0.7: '#62B4E7', 1.:'#4584F0'})
    plotitem.pcolor_cmap = white_green_cmap
    #plotitem.add_colorbar = True
    plotitem.plot_var = Pressure  # defined above
    #plotitem.plotstyle = '-o'
    #plotitem.color = 'r'
    # For AMR patches and cell edges (# REMEMBER TO CHANGE amr_contour_show TOO)
    plotitem.amr_patchedges_show = [0,0,0,1]  #[0,0,0,0,1] #[0,0,0,0,0,1]
    plotitem.amr_celledges_show = [1,1,1,0] #[1,1,0,0,0] #[1,1,1,1,0,0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True 
    
    # Add contours as well
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = Pressure
    plotitem.contour_levels = np.linspace(90000,230000,30) 
    #plotitem.contour_nlevels = 10
    #plotitem.contour_min = 91000.0
    #plotitem.contour_max = 290000.0
    #plotitem.amr_patchedges_show = [0,0,1]
    #plotitem.amr_celledges_show = [1,1,0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True 
    plotitem.amr_contour_colors = ['b','#3C3C3C','k']
    plotitem.amr_contour_show = [0, 0, 0, 1]

    
    plotaxes.afteraxes = aa
    
    # Figure for Pressure (Schlieren)
    plotfigure = plotdata.new_plotfigure(name='Pressure schlieren', figno=9)
    plotfigure.kwargs = {'figsize':[8,3.7], 'tight_layout':True}
    #plotfigure.kwargs = {'figsize':[8,8], 'tight_layout':True} # For colorbar output

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('Pressure')
    plotaxes.xlimits = [-0.04,0.04] 
    plotaxes.ylimits = [0.001,0.035]
    plotaxes.title = 'Pressure sclieren'
    plotaxes.scaled = True      # so aspect ratio is 1
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.schlieren_cmin = 500 #2000 #500 #20
    plotitem.schlieren_cmax = 30000 #3500 #25000 #30000
    plotitem.add_colorbar = True
    plotitem.plot_var = Pressure  # defined above
    # For AMR
    plotitem.amr_patchedges_show = [0,0,0,0,1]
    plotitem.amr_celledges_show = [0,0,0,0,0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True 
    
    plotaxes.afteraxes = aa

    
    # Figure for Pressure slice
    # -------------------
    
    plotfigure = plotdata.new_plotfigure(name='Pressure slice', figno=6)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    # Axes for m vs Pa or cm vs PSI
    #plotaxes.xlimits = [-0.03,0.03] #[-3,3] #[-8.5,16] #'auto' -16
    #plotaxes.ylimits = [0.00000,300000]
    plotaxes.xlimits = [-3.0,3.0]
    plotaxes.ylimits = [-20,30]
    plotaxes.title = 'Pressure slice'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        aux = current_data.aux

        ij = find((y <= dy/2.) & (y > -dy/2.))
        x_slice = ravel(x)[ij]
        gamma_slice = ravel(aux[0,:,:])[ij]
        pinf_slice = ravel(aux[1,:,:])[ij]
        rho_slice = ravel(q[0,:,:])[ij]
        momx_slice = ravel(q[1,:,:])[ij]
        momy_slice = ravel(q[2,:,:])[ij]
        ene_slice = ravel(q[3,:,:])[ij]
        P_slice = (gamma_slice - 1.0)*(ene_slice - 0.5*(momx_slice**2 + momy_slice**2)/rho_slice)
        P_slice = P_slice - gamma_slice*pinf_slice
        # Convert to Psi and centimeters
        P_slice = P_slice*0.000145038 - 14.6959488
        x_slice = 100*x_slice
        return x_slice, P_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = '-kx'
    plotitem.kwargs = {'markersize':3}
    
    plotaxes.afteraxes = aa1DPSIcm
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                      # print figures
    plotdata.print_format = 'png'                  # file format
    plotdata.print_framenos = [32,57,68,74,88,108] # list of frames to print 'all' for all frames
    plotdata.print_fignos = [4] #'all'             # list of figures to print
    plotdata.html = True                           # create html files of plots?
    plotdata.html_homelink = '../README.html'      # pointer for top of index
    plotdata.latex = True                          # create latex file of plots?
    plotdata.latex_figsperline = 2                 # layout of plots
    plotdata.latex_framesperline = 1               # layout of plots
    plotdata.latex_makepdf = False                 # also run pdflatex?

    return plotdata

    
