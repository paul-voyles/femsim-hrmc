If you use FEMSIM or HRMC, please cite this [publication](http://www.sciencedirect.com/science/article/pii/S001046551630385X):

    J. J. Maldonis, J. Hwang, P. M. Voyles, FEMSIM + HRMC: Simulation of and structural refinement using fluctuation electron microscopy for amorphous materials. Comput. Phys. Commun. (2016), doi:10.1016/j.cpc.2016.12.006.

This publication can also be used for further documentation.

## Input and Output Parameters and Files

###  Compiling
All code is written in FORTAN and can be compiled using an Intel Fortran 90 compiler that supports MPI (for example, Intel’s mpif90).  Compilation has only been tested on linux.  Compilation should be done using GNU Make and the makefile provided in the `src` directory.  A make target is required, and can be either `hrmc` or `femsim`.  An optional keyword `debug` will compile the code with debugging enabled.

For example, `make hrmc debug –C src` (run from the root directory) will create an executable with the name `hrmc`, compiled without optimizations and with debugging enabled.

Once compiled, the code can be run with `./executable basename paramfile.in` where `executable` is either `hrmc` or `femsim`, `basename` specifies a keyword that is used as part of the output filenames to prevent output files from being overwritten when multiple simulations are run in the same directory, and `paramfile.in` is the name of a parameter file.

###  Input files

The code is contained in the following files:
* `eam.f90` and `eam_fs.f90` define the EAM/alloy and EAM/Finnis-Sinclair potential energy functionality, respectively
  * The user will need to specify which functional form should be used by including either `eam.f90` or `eam_fs.f90` within the source files in the makefile.
* `fem.f90` contains the code relevant to the FEMSIM calculation, including its usage within HRMC
* `globals.f90` defines a few global variables (e.g. the MPI variables) for the entirety of the code
* `model.f90` defines the code relevant to handling the atomic model
* `read_inputs.f90` parses the parameter file and sets relevant variables
* `rmc_functions.f90` defines functions necessary for the HRMC looping schema
* `hrmc.f90` contains the core code that runs HRMC or FEMSIM, depending on the make command
* `scattering_factors.f90` defines the electron scattering factors used within the FEMSIM calculation

Additional input files for FEMSIM include:
* an atomic model file in XYZ format
  * There is one additional formatting requirement:  The comment line must start with three floats (separated by a space) indicating the size in Angstroms of the atomic model.  The atomic model is required to be cubic.
* a file containing a list of `k` points for which `V_sim(k)` will be evaluated
  * The first row is comment line, followed by any number of `k` values, each on a separate line.
* a submit file for submitting the code to a cluster (optional)

Additional input files for HRMC include:
* an atomic model file in XYZ format
  * There is one additional formatting requirement:  The comment line must start with three floats (separated by a space) indicating the size in Angstroms of the atomic model.  The atomic model is required to be cubic.
* an experimental data file containing `V(k)` information
  * The first row is a comment line, followed by any number of rows of data points.
  * The data is space-delimited and must be three columns long:  1) `k`, at which `V_sim(k)` will be evaluated; 2) `V_exp(k)` for comparison to `V_sim(k)`; and 3) `sigma(k)^2`.
* an EAM/alloy or EAM/Finnis-Sinclair potential file in standard format, with one float per line
  * the helper file `utils/reformat_potential.py` may be helpful for putting a single float on each line
* a submit file for submitting the code to a cluster (optional)
  * An example using the SLURM scheduling system is included.

The parameter file must have the following format on each line:

1. comment line, ignored by the program

2. input atomic model filename

3. input experimental data filename

4. `Q = 0.61/R` (defined by the Rayleigh resolution criterion)

5. number of rotations: `nphi`, `npsi`, `ntheta`

  a. `phi`, `psi`, and `theta` follow the Goldstein convention, and `nphi`, `npsi`, and `ntheta` represent the number of angles about their respective axes.

  b. The angles are chosen such that the unit sphere is uniformly sampled.

  c. `ntheta` is the number of angles about the z-axis and should be set to `1` for all practical purposes.

  d. Note that `phi` ranges from `0 - 2pi`, `psi` ranges from `0 - 2pi`, and `theta` ranges from `0 - pi`.

  e. In the FEMSIM example, `nphi = 1`, `npsi = 40`, `ntheta = 20`.

The following lines can be omitted for FEMSIM but are required for HRMC (continued on newlines after the FEMSIM parameters above):

6. thickness scaling factor: `beta = 3 * experimental_sample_thickness / simulation_side_length`

7. input potential filename

8. `starting_step`, `ending_step`

  a. `starting_step` is the step number at which the simulation will start and is useful for continuations of a previous simulation.

  b. `ending_step` is the step number at which the simulation will terminate and is most useful when a cluster’s queuing system has a maximum job duration.

9. `temperature`, `max_move`, `decrement`
  a. `temperature` is the temperature of the simulation at `starting_step = 0`.  

  b. `max_move` is the maximum distance an atom is allowed to move at `starting_step = 0`.

  c. `decrement` is the number of steps that are run before temperature and max_move are decremented via an inverse power law.

10. random number generator seed (integer only)

11. weighting factor, `alpha` , between `E` and `chi`<sup>2</sup>

12. number of atom species in the atomic model

13. The next N lines are a list of N hard sphere cutoff distances, where N is the number of atom species. The column/row matrix has the same order (from left to right and top to bottom) as the atom species in the potential file, and values on each line must be separated by a space.  When an atom is randomly moved in the HRMC algorithm, the move will be immediately retried without calculation of the energy or `chi`<sup>2</sup> if moved atom is within the hard sphere cutoff distance of another atom.  We suggest using the minimum interatomic distance from `g(r)` data for the hard sphere cutoff values.

###  Output files

The `basename` keyword in the following names is a placeholder for a command line input that prevents output files from being overwritten when multiple simulations are run in the same directory.  In the case that HRMC is submitted to a cluster, `basename` can be replaced by the job ID.

* `stdout` displays step-by-step information, including whether an atom move was accepted or rejected and the values that pertain to the acceptance criterion.
  * Written to every step.
* `stderr` displays all error related information.
  * Written to when an error occurs.
* `acceptance_rate_basename.txt` contains a single comment line followed by columns of step numbers, in increments of 100, and acceptance rates averaged over the previous 100 steps.
  * Written to every 100 steps.
* `chi_squared_basename.txt` contains a comment line followed by columns of step number, `chi`<sup>2</sup>, and `E`.
  * Written to every time an atom move is accepted.
* `model_update_basename_N.xyz` contains the model file in XYZ format at step N.
  * Written once every 1000 steps.
* `model_final_basename.xyz` contains the model file in XYZ format at the end of the simulation.
  * Written once at the end of the simulation.
* `vk_initial_basename.txt` contains the simulated `V(k)` for the starting model.
  * Written once at the beginning of the simulation.
* `vk_final_basename.txt` contains the simulated `V(k)` for the final model.
  * Written once at the end of the simulation.
