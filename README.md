# Group Project: Solving the Two Dimensional Poisson Equation

## Introduction

The program solves Poisson's equation, <img src="/tex/a728fdceeb7edb6d08c57ab25335a6f7.svg?invert_in_darkmode&sanitize=true" align=middle width=133.56565365pt height=26.76175259999998pt/> in two dimensions through either a parallelization or serial verison of the code. For the former, we accomplish parallelization via OpenMP.

## Theory
























## Methods (Methods are all copied from Assignment pard of README file)

In this project you will solve Poisson's equation, <img src="/tex/a728fdceeb7edb6d08c57ab25335a6f7.svg?invert_in_darkmode&sanitize=true" align=middle width=133.56565365pt height=26.76175259999998pt/> ,where <img src="/tex/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999991pt height=22.831056599999986pt/> is the electrostatic potential and <img src="/tex/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode&sanitize=true" align=middle width=8.49888434999999pt height=14.15524440000002pt/>
is the source charge density. Towards the end of the course we will learn how
to solve this via a "relaxation" method, but here we will transform the
differential equation into a matrix problem in a given basis. The main
numerical work will be Fourier integrals.

We will work in two dimensions, so that the full equation becomes
<p align="center"><img src="/tex/31659c8b0dc6dfd77f3ad90c50fb13cf.svg?invert_in_darkmode&sanitize=true" align=middle width=254.44296074999997pt height=40.11819404999999pt/></p>

We will place this in a metal rectangular box with sides of length <img src="/tex/9dbcef13f3e6981dfe63f653112a933f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.64161584999999pt height=22.465723500000017pt/> and
<img src="/tex/ac6244c7cc0673f3a8d51603f75bcffe.svg?invert_in_darkmode&sanitize=true" align=middle width=18.26684969999999pt height=22.465723500000017pt/>. This means the gradient of the potential at the surfaces must be zero,
that is, if we place the walls at <img src="/tex/8436d02a042a1eec745015a5801fc1a0.svg?invert_in_darkmode&sanitize=true" align=middle width=39.53182859999999pt height=21.18721440000001pt/> and <img src="/tex/9dbcef13f3e6981dfe63f653112a933f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.64161584999999pt height=22.465723500000017pt/>, and at <img src="/tex/a42b1c71ca6ab3bfc0e416ac9b587993.svg?invert_in_darkmode&sanitize=true" align=middle width=38.78604674999999pt height=21.18721440000001pt/> and <img src="/tex/ac6244c7cc0673f3a8d51603f75bcffe.svg?invert_in_darkmode&sanitize=true" align=middle width=18.26684969999999pt height=22.465723500000017pt/>, then
<p align="center"><img src="/tex/121b2aba0948920e666965710a533721.svg?invert_in_darkmode&sanitize=true" align=middle width=164.84245965pt height=92.0098938pt/></p>

This is easily accomplished by expanding the potential in basis functions
which automatically satisfy these conditions, that is cosines:
<p align="center"><img src="/tex/446111e69ab8e69beb7d8ffa09e7eff6.svg?invert_in_darkmode&sanitize=true" align=middle width=368.4128217pt height=49.73538075pt/></p>
These basis functions have been chosen so as to be already normalized.

In this basis, the problem becomes almost trivial:
<p align="center"><img src="/tex/132295f99071c92887c636369deed287.svg?invert_in_darkmode&sanitize=true" align=middle width=261.1920333pt height=49.315569599999996pt/></p>

to be solved for <img src="/tex/3343d1e4776a8c1cc78c3aa3e1c3c557.svg?invert_in_darkmode&sanitize=true" align=middle width=30.808809899999993pt height=14.15524440000002pt/> and where the source term is now

<p align="center"><img src="/tex/4ec3243c4b95787028f2d53b801683ad.svg?invert_in_darkmode&sanitize=true" align=middle width=447.1595766pt height=44.749102199999996pt/></p>

This two dimensional integral is the main computational work. For this you can
***and should*** use your `quadrature` module from assignment two. The number
of sample points for the quadrature ***does not*** need to be the same as the
number of samples for the Fourier coefficients. You may hard code the number
of quadrature sample points. If your quadrature subroutines are general enough
you just need to copy the `quadrature.f90` file into this repository and you
wont even need to make changes to the source code. However you will need to
modify the `makefile` in this repository to include the correct dependencies.

For the source term, we will use the function
<p align="center"><img src="/tex/59c9ff9b238c60d6db45edb314c713c8.svg?invert_in_darkmode&sanitize=true" align=middle width=343.05987539999995pt height=42.07871745pt/></p>

Your code should receive, via a `namelist` file given as an argument, the
following:

* The dimensions of the box, <img src="/tex/9dbcef13f3e6981dfe63f653112a933f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.64161584999999pt height=22.465723500000017pt/> and <img src="/tex/ac6244c7cc0673f3a8d51603f75bcffe.svg?invert_in_darkmode&sanitize=true" align=middle width=18.26684969999999pt height=22.465723500000017pt/> (which may be different)
* The total charge <img src="/tex/3ee3f93b7a51e5719c84fadc68137817.svg?invert_in_darkmode&sanitize=true" align=middle width=15.05143034999999pt height=14.15524440000002pt/>.
* The center of the charge distribution <img src="/tex/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> and <img src="/tex/14adeddbb1889c9aba973ba30e7bce77.svg?invert_in_darkmode&sanitize=true" align=middle width=14.61197759999999pt height=14.15524440000002pt/>
* The widths <img src="/tex/baed07e6cbaba3e37ce167d64db1675d.svg?invert_in_darkmode&sanitize=true" align=middle width=19.935847799999987pt height=22.465723500000017pt/> and <img src="/tex/235409eac7aaaacd0cc98f7e38ea1ecd.svg?invert_in_darkmode&sanitize=true" align=middle width=19.561081649999988pt height=22.465723500000017pt/>
* The max index <img src="/tex/58a0abdbe1ebe952d7251c442f2dee8d.svg?invert_in_darkmode&sanitize=true" align=middle width=37.46589329999998pt height=22.465723500000017pt/>, which is the maximum value for both indices <img src="/tex/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode&sanitize=true" align=middle width=14.433101099999991pt height=14.15524440000002pt/> and <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/>.

You may add additional `namelists` to the input file (e.g. name of the output
file, number of sample points for the quadrature, step size in each direction
when writing results for plots), just be sure to give default values to all
variables in all `namelists`.

Create a jupyter notebook to show your resulting potential <img src="/tex/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999991pt height=22.831056599999986pt/> in a two
dimensional plot. Either through iso-potential curves (`plt.contour()`) or by a
color plot (`plt.imshow()`). I'll leave up to you to decide how to better
present your results and to look at the necessary documentation.

### Serial

# Integration Algorithm
We selected the Boole's algorithm 




### Parallelization

To utilize the efficiency of OpenMp, the program puts the parallel codes in the data file generating part of `read_write.f90`, them the task of generating data file is available to break into n (determined by how many threshold) parallel processes to save time.

## Results




### Convergence criteria: how to know when you have a reasonable n_max.

The n_max can be chosen in different ways depending on what you are looking for:
If your goal is just to get an understanding of what the structure of the potential in the box is, an `n_max` between 20 - 40 would suffice.
In the case you are looking for precision, an n_max no greater than 70 should be useful in getting a good approximation.





## Conclusion

### What to expect from user.

* The program contributes a simpiler way to show the process of solving Poisson's equation with codes instead of equations or methods, which allows the user to understand where and how those values are used and contributed to the results. In the meantime, the program is a convenient tool for analyzing the difference from 2D and 3D potential graphs by slightly switching the initial conditions in the `.namelist` file.

* 2D and 3D graphs in `.ipynb` are able to exhibit a clear view of distribution of potential.

* Be familiar with the power ot OpenMP that the efficiency and time the program takes for computation by changing the number of threshold.

### How should the program work

The program relies on fortran language to code the formula and system required to run the program which returns a output of database including positions and potentials into the code in jupyter notebook which is already prepared to draw the graph we expect.



### Program Contents

The program includes:

* `box_parameters.namelist` file for convenient to change initial conditions.
* `graoup_project_analysis.ipynb` shows the results and graphs in jupyter notebook.
* `main.f90`, the overview of program with calling all necessary subroutines and functions from other files in order to generating data files for graphing.
* `makefile`, a shortcut code for executing the program.
* `potential.f90` is responsible solving Poisson's equation for a charge distribution held within a 2D conducting box. Ultimately, this file finds the electrostatic potential as a function of x and y.
* `quadrature.f90` includes 2 functions used to evaluate a volume integral with the methods of Booles' quadrature.
* `read_write.f90` is the file used to read initial values for computing the potential and generate data files for graphing in the '.ipynb' file.

### Instructions for Usage

Determine the initial condition inside `.namelist` file (or the program will take initial values from the `read_write.f90` file). Due to the integration algorithm of Boole's rule, the value of `n_max` has to be 4n+1. And 0 value in `center` or `width` would cause invalid output values.

Once the initial conditions are setting, open the terminal and go to the direction where the program is, write `make` and `./poisson box_parameters.namelist` (or `./poisson` to use default initial condition) in order to compiling and executing the program. In generally it will takes around 15 seconds of executing and the time will be showd on the screen of terminal.

When the `.dat` file is generated, the jupyter notebook file `.ipynb` is available to graph the distribution of potential with both 2D and 3D graphs. And the jupyter notebook requires the `numpy` and `matplotlib` library for graphing.

Note: write `make clear` in the terminal is available to clear all unecessary files which used to generating the execution file.


