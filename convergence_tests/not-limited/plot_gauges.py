
from pylab import *
from clawpack.visclaw.data import ClawPlotData
from scipy import interpolate

# Aux Parameters
gammawat = 7.15
pinfwat = 300000000.0
rhow = 1000.0
gaugeno = 1 #2

plotdata = ClawPlotData()

# Folder from out or out2 differ in that the Riemann solver used uses the Lagrangian transformation
# on all of water or only in the interface respectively

# CHOOSE out or outmap in the outdir to choose from non-mapped or mapped version of code
outstr = '_out'

#plotdata.outdir = '../_output'   # set to the proper output directory
#g = plotdata.getgauge(gaugeno)
#p = zeros(g.t.size)
#p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
#p = 0.001*p # Convert to KPa
#tt = g.t*1000000 # Convert to microsec
#plot(tt, p, '-g', label="Level New; ", linewidth=2)

plotdata.outdir = outstr +'_conv_40x20_lvl6_refrat2-2-2-2-2'   # set to the proper output directory
g = plotdata.getgauge(gaugeno)
p = zeros(g.t.size)
p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
p = 0.001*p # Convert to KPa
tt = g.t*1000000 # Convert to microsec
plot(tt, p, '-k', label="Level 6", linewidth=1)

plotdata.outdir = outstr +'_conv_40x20_lvl5_refrat2-2-2-2'   # set to the proper output directory
g = plotdata.getgauge(gaugeno)
p = zeros(g.t.size)
p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
p = 0.001*p # Convert to KPa
tt = g.t*1000000 # Convert to microsec
plot(tt, p, '-b', label="Level 5", linewidth=1)

# Load and plot ngauge from 40x20 grid AMR run with level 3 AMR and refrinement ratio [2,2]
plotdata.outdir = outstr +'_conv_40x20_lvl4_refrat2-2-2'   # set to the proper output directory
g = plotdata.getgauge(gaugeno)
p = zeros(g.t.size)
p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
p = 0.001*p # Convert to KPa
tt = g.t*1000000 # Convert to microsec
plot(tt, p, '-r', label="Level 4", linewidth=1)

# Load and plot ngauge from 40x20 grid AMR run with level 3 AMR and refrinement ratio [2,2]
plotdata.outdir = outstr +'_conv_40x20_lvl3_refrat2-2'  # set to the proper output directory
g = plotdata.getgauge(gaugeno)
p = zeros(g.t.size)
p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
p = 0.001*p # Convert to KPa
tt = g.t*1000000 # Convert to microsec
plot(tt, p, '-g', label="Level 3",linewidth = 1)

## Load and plot ngauge from 40x20 grid AMR run with level 3 AMR and refrinement ratio [2,2]
plotdata.outdir = outstr +'_conv_40x20_lvl2_refrat2-2'   # set to the proper output directory
g = plotdata.getgauge(gaugeno)
p = zeros(g.t.size)
p = (gammawat - 1.0)*(g.q[3,:] - 0.5*(g.q[1,:]*g.q[1,:] + g.q[2,:]*g.q[2,:])/g.q[0,:]) - gammawat*pinfwat
p = 0.001*p # Convert to KPa
tt = g.t*1000000 # Convert to microsec
plot(tt, p, '-y', label="Level 2", linewidth=1)

plt.xlabel('time $(\mu s)$', fontsize=16)
plt.ylabel('pressure $(KPa)$', fontsize=16)
plt.axis([60,100,90,300])
plt.legend(loc = 'upper left')
plt.show()


#

#title("Gauges on transducer (red near center, black near edge)")


#figure(figsize=(12,6))
#plot(g0.t, p_average, 'b')
#plot(tt, pt, 'r')
#title("Average pressure on transducer")

# Gauges locations from setrun
    #gauges.append([0, -0.01, 0, 0., 1e9])
    #gauges.append([1, -0.01, 0.005, 0., 1e9])
    #gauges.append([2, -0.01, 0.01, 0., 1e9])
    #gauges.append([3, 0.0, 0, 0., 1e9])
    #gauges.append([4, 0.0, 0.005, 0., 1e9])
    #gauges.append([5, 0.0, 0.01, 0., 1e9])
    #gauges.append([6, 0.01, 0, 0., 1e9])
    #gauges.append([7, 0.01, 0.005, 0., 1e9])
    #gauges.append([8, 0.01, 0.01, 0., 1e9])

