# [Recursive Weighted Myriad Based Filters and Their Optimizations](http://ieeexplore.ieee.org/document/7457715/)

[Juan Marcos Ramirez](juanra@ula.ve) and [Jose Luis Paredes](https://www.eecis.udel.edu/~paredesj/)

IEEE Transactions on Signal Processing

## Abstract
This paper proposes two new recursive filtering structures based on the nonlinear myriad operator. First, we develop the general class of *recursive weighted myriad* (RWMy) filters as a robust filtering structure against impulsive noise that includes, as particular cases, a normalized version of the linear infinite impulse response (IIR) filter, the recursive mode-type filter, and the non-recursive weighted myriad filter. Secondly, considering the fact that the additive noise that contaminates the previous filter's outputs is no longer impulsive, we introduce a novel class of *recursive hybrid myriad* (RHMy) filters whose structure gathers the advantages of both the weighted myriad and the weighted mean in a single cost function to be minimized. Some properties of the RHMy filter are derived and a fast algorithm to compute the RHMy filter output is proposed. Furthermore, adaptive algorithms for designing the proposed recursive structures, under the equation error formulation framework, are developed. Finally, extensive numerical simulations are shown to evaluate both the iterative update of the adaptive algorithms and the performance of the proposed recursive filters against impulsive noise.

## Suplementary Materials

### How to run the code

Download and uncompress the `RecursiveMyriadBasedFilters` folder. To generate Figures and Tables in the paper, under **MATLAB** environment, navigate to the `RecursiveMyriadBasedFilters` folder and follow the instructions described below

#### Figure 4

To generate Figures 4(a)-4(e) run, in MATLAB, 

	>> FigTrainingBehavior.m

To generate Figures 4(f) run 

	>> FigComputationalCost.m

#### Figure 5

	>> FigCleanChirpFiltering.m

#### Figure 6

	>> FigNoisyChirpFiltering.m

#### Table I

	>> TableNoisyChirp.m

#### Table II

	>> TableImageDenoising.m

#### Figure 8

	>> FigPLCSignals.m 

#### Figure 9

	>> FigPLCMetrics.m

### Platform

* Scientific Fedora 21 Operating System, MATLAB R2013a


## License

This code package is licensed under the GNU GENERAL PUBLIC LICENSE (version 3) - see the [LICENSE](LICENSE) file for details


## Author

Juan Marcos Ramírez Rondón. Escuela de Eléctrica. Facultad de Ingeniería. [Universidad de Los Andes](http://www.ula.ve). Mérida, 5101, Venezuela. 


### Contact

[Juan Marcos Ramirez](juanra@ula.ve)

## Date

September 14, 2016
