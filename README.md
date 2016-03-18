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

* **convergence_tests:** contains convergence tests by comparing the output data at given gauges for many different simulations using different levels of AMR refinement. The internal folders: not-limited and limited, contain the convergence tests before and after applying the limiters developed in the paper. Each of the folders contain subfolders with the gauge output data for several simulations at different levels of AMR refinement. Both folders also contain a python script to produce the convergence test plots from the paper. 

### Changing the output
Although changing the output might require getting more involved with the code, there are some simple tweaks that will allow the user to see different output. The main code files to edit are setrun.py and setplot.py

**setrun.py**

* Change domain boundaries:

```
    clawdata.lower[0] = -0.03                 # xlower
    clawdata.upper[0] = 0.03                  # xupper
    clawdata.lower[1] = 0.000000e+00          # ylower
    clawdata.upper[1] = 0.020000e+00          # yupper
```

* Change grid size (from 600x60 to new values):

```
    clawdata.num_cells[0] = 600
    clawdata.num_cells[1] = 60
```

* Change output times:

```
    clawdata.num_output_times = 500 # number of output frames
    clawdata.tfinal = 0.0002 # final time in seconds
```

**setplot.py**

* Change main properties of figure 7 (the one showing in the paper):

```
    # Pressure contour(2D) and pressure slice(1D) in one figure
    plotfigure = plotdata.new_plotfigure(name='Contour & Slice', figno=7)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes() 
    plotaxes.xlimits = [-3,3] 
    plotaxes.ylimits = [-20,30]
    plotaxes.title = 'Pressure'    
    plotaxes.afteraxes = MirrorPressurecontour_N_Pressureslice    
 ``` 
 
* Change plotting output:

```
    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos =  [75,150,158,174,212,336] #list of frames to print. Use 'all' to print all #
    plotdata.print_fignos = 'all'            # list of figures to print
