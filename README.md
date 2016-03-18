## Numerical methods for interface coupling of compressible and almost incompressible media.

This is the code to supplement the paper: "Numerical methods for interface coupling of compressible and almost incompressible media.". It is an extended and more general version of the methods developed in [BBB_experiment](https://github.com/maojrs/BBB_experiment). It employs 2D axisymmetric Euler equations coupled with a Tammann equation of state, and it uses Lagrangian/Eulerian coupling along interfaces between compressible and almost incompressible media. In addition, it implements mapped grids to model more complex interface geometry, and it uses adaptive mesh refinement (AMR) to provide better accuracy while diminishing computational cost.

The code available here serves as an implemented example of the methods described in the paper, and it reproduces its main figures. It allows the user to easily modify the code, so it can be employed in other applications by research scientists. 

The code is designed to work with the Clawpack 5 software. For more information visit [the Clawpack webpage](http://www.clawpack.org/ ). The Clawpack software is open source, and the code is openly available at [Clawpack on GitHub](https://github.com/clawpack/clawpack). Some of the main algorithms used in this project can be found in the book [Finite Volume Methods for Hyperbolic Problems](http://depts.washington.edu/clawpack/book.html).

### Software dependencies
**All the required software is open source.**
* gfortran 4.8.1 
* python 2.7.5
* clawpack 5.3.1 : [Clawpack](http://www.clawpack.org/ )
* git: to clone this repository

--Clawpack 5.0 or higher is required, other versions of gfortran and python might work as well.

**Operating system.**
* Linux
* Mac OS X

### How to run this code?
1. **Install Clawpack 5.3.1:**
  * Go to: http://www.clawpack.org/installing.html#installation-instructions
  * Follow the download and installation intructions on the section: Install all Clawpack packages and Set environment variables.
  * Test your installation running an example for Classic Clawpack on the section: Testing your installation.
  * If your installation works, you already have python and gfortran installed.

2. **Clone this repository to your local machine:**

 ```
    git clone https://github.com/maojrs/Interface_Euler_AMR.git
 ```

3. **Run the code**
Go to the folder: *code_cartesian* or *code_mapped* and run in a terminal:

 ```
    make .plots
 ```

The code will produce two new folders: _output and _plots. The first one contains all the output files, while the latter one contains the plots and interactive visualization apps.

### Folder organization

* **code_cartesian:** contains the code to run the 2D axisymmetric simulation of Euler equations with AMR on a Cartesian grid with a rectangular interface between air and water. The internal folders: _initfiles, _plots_paper and rp, contain the initial condition data for the incoming shock wave, the plots produced included in the paper, and the code for the corresponding Riemann solvers, respectively.

* **code_mapped:**  contains the code to run the 2D axisymmetric simulation of Euler equations with AMR on a mapped grid with a circular interface between air and water. The mapped grid can be modified by the user to implement different mappings. The mapping included in mapc2p.py and mapc2p.f90 not only allows a circular inclusion but also a ring inclusion. The internal folders: _initfiles, _plots_paper and rp, contain the initial condition data for the incoming shock wave, the plots produced included in the paper, and the code for the Riemann solvers in the mapped grid, respectively.

* **convergence_tests:** contains convergence tests by comparing the output data at given gauges for many different simulations using different levels of AMR refinement. The internal folders: default-limiter and new-limiter, contain the convergence tests before and after applying the limiters developed in the paper. Each of the folders contain subfolders with the gauge output data for several simulations at different levels of AMR refinement. Both folders also contain a python script to produce the convergence test plots from the paper. 

### Changing the output
Although changing the output might require getting more involved with the code, there are some simple tweaks that will allow the user to see different output. The main code files to edit are setrun.py and setplot.py. These files can be found in the cartesian and in the mapped grid code.

**setrun.py**

* Change computational domain boundaries:

```
    clawdata.lower[0] = -0.05          # xlower
    clawdata.upper[0] = 0.05           # xupper
    clawdata.lower[1] = 0.000000e+     # ylower
    clawdata.upper[1] = 0.05000e+00    # yupper
```

* Change initial grid size (from 40x20 to new values):

```
    clawdata.num_cells[0] = 40
    clawdata.num_cells[1] = 20
```

* Change output times:

```
    clawdata.num_output_times = 150 # number of output frames
    clawdata.tfinal = 0.0002        # final time in seconds
```

* Add or remove output gauges at position x, y.
```
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    gauges.append([0, -0.01, 0, 0., 1e9])
    gauges.append([1, -0.01, 0.005, 0., 1e9])
    gauges.append([2, -0.01, 0.01, 0., 1e9])
    gauges.append([3, 0.0, 0, 0., 1e9])
    gauges.append([4, 0.0, 0.005, 0., 1e9])
    gauges.append([5, 0.0, 0.01, 0., 1e9])
    gauges.append([6, 0.01, 0, 0., 1e9])
    gauges.append([7, 0.01, 0.005, 0., 1e9])
    gauges.append([8, 0.01, 0.01, 0., 1e9])
```

* Change AMR refinement level and refinement ratios
```
    amrdata.amr_levels_max = 4
    # List of refinement ratios at each level (length at least amr_level_max-1)
    # Note changing refinement ratio might affect stability at interface
    amrdata.refinement_ratios_x = [2, 2, 2]
    amrdata.refinement_ratios_y = [2, 2, 2]
    amrdata.refinement_ratios_t = [2, 2, 2]
 ```
 
* Force AMR refinement regions
```
    # to specify regions of refinement append lines of the form
    # regions.append([minlevel,maxlevel,t1,t2,x1,x2,y1,y2])
    regions.append([4,4,0,1e9,-0.0155,0.0155, 0.0, 0.0155])
```


**setplot.py**

* Change main properties of figure 4 (the one showing in the paper):
```
    plotaxes = plotfigure.new_plotaxes('Pressure')
    plotaxes.xlimits = [-0.04,0.04] 
    plotaxes.ylimits = [0.001,0.035]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1   
 ``` 
 
* Change colormap and range of the figure
```
plotitem.pcolor_cmin = 90000
    plotitem.pcolor_cmax = 230000
    #plotitem.pcolor_cmap = colormaps.white_blue
    white_green_cmap = colormaps.make_colormap({0.:'w', 0.35: '#AAFFEF', 0.7: '#62B4E7', 1.:'#4584F0'})
```
 
* Change AMR patches and grids shown in figure
```
    # Each position of the array corresponds to an AMR refinement level.
    plotitem.amr_patchedges_show = [0,0,0,1] 
    plotitem.amr_celledges_show = [1,1,1,0] 
```
 
* Change contours shown on top of the pcolor plot
```
    plotitem.contour_levels = np.linspace(90000,230000,30)
    ...
    # Show contours only for highest level 
    plotitem.amr_contour_show = [0, 0, 0, 1]
```
 
* Change plotting output:
```
    plotdata.printfigs = True                      # print figures
    plotdata.print_format = 'png'                  # file format
    plotdata.print_framenos = [32,57,68,74,88,108] # list of frames to print 'all' for all frames
    plotdata.print_fignos = [4]                    # list of figures to print 'all' for all figures
```
