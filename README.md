# Parameter file input
Lines below start with the line number in the parameter file. See param_file.in for an example. All necessary values will be encased in <> below with a description of each further below.

1 A comment line. It doesn't matter what, just make sure it's there.

2 \<modelfile>

3 \<v(k) file>

4 \<eam potential>

5 \<starting step #> \<ending step #>

6 \<starting temperature> \<starting maxmove> \<how often to decrement temp and max move>

7 \<number of elements in the model>

8-... \<hard sphere cutoff distances>

9 \<alpha>

10 \<thickness scaling factor>

11 \<rotations, three values>

12 \<Q>

\<modelfile>:  The starting modelfile for the HRMC simulation. The format is not the standard xyz file format, see examples included in repo.

\<v(k) file>:  The experimental V(k) data to fit against. See examples included in repo for formatting.

\<eam potential>:  An EAM potential file for energy calculations.

\<starting step #>:  The starting step # to be output in the chi_squared.txt file. For easy of graphing.

\<ending step #>:  The job will run until this step number is hit.

\<starting temperature>:  The starting temperature, assuming the step # is 0. Do not modify this if you are continuing a simulation -- just make sure your starting step # is correct.

\< starting maxmove>:  The starting maximum move distance an atom is allowed to move, , assuming the step # is 0. Do not modify this if you are continuing a simulation -- just make sure your starting step # is correct.

\<how often to decrement temp and max move>:  The number of steps between decreasing the temperature and maximum move distance for simulated annealing. I set this to approximately 1/20 of the final simulation # of steps. You can calculate a cooling rate with respect to this and gauge whether you are cooling too quickly or not. You want to be able to explore as much of the potential energy landscape as you can bear (time is always an issue).

\<number of elements in the model>:  Self explanitory.

\<hard sphere cutoff distances>:  This will be a series of (# of elements) x (# of elements) lines that contains the hard-sphere cutoff distances you want to use. No atom move is allowed that puts atoms of specific types closer together than this. Different values can be input for specific atomic species pairs. The order from left to right and top to bottom should be printed shortly after the simulation begins (that should be correct).

\<alpha>:  See description below. This deserves its own heading.

\<thickness scaling factor>:  The ratio of the experimental sample thickness to the model thickness, multiplied by 3. The 3 accounts for approximations in the experimental data simulation algorithm.

\<rotations, three values>:  Three values for determining how to rotate the primary model for calculations of V(k). Use 1 40 20, or read more below.

\<Q>:  Q = 0.61/(FEM experimental resolution). Make sure the units match those in the modelfile. Use Angstroms.

### Other notable parameters not in the parameter file
There is a "do while" loop in hrmc.f90 that specifies how long to run for. I recommend running indefinitely, and manually cancelling the job when you decide it is done based on the contents of chi_squared.txt. However, if your cluster only lets you run for a certain amount of time, you may need to modify this.

# How the simulation works
HRMC stands for Hybrid Reverse Monte Carlo. Monte Carlo (MC) is the simulation method. Reverse means we refine against experimental data. Hybrid means we refine against an energetic potential as well as experimental data.

MC is a statistical method. Ideally, we want to run multiple HRMC simulations and analyze statistically. However, time is an issue. It can be helpful to encourage the results to be as "different" as possible. One method of doing this is starting with a crystalline and an amorphous structure as the initial model. How you create these starting models is up to you.

The simulation of the experimental data, in our case the variance of a set of intensity vs k curves, is detailed in a forthcoming paper.

### Rotations
The number of rotations in the model can vary, we typically set it to 1 40 20 = 211 rotations. That means the intensity calculation is being done that many times, once for each rotation, and then combined into the V(k) curve for experimental comparison.

Note that the first rotation should be 0 0 0, so that the first model in the rotation array (mrot) is identical to the model in hrmc.f90.

### Pixels
The model must obey periodic boundary conditions, and therefore must be rectangular. It is futhermore necessary that the model is cubic, not just rectangular. Every atom in the model must contribute once and only once to the intensity calculation (per rotation). With round pixels (which correctly represents the exeriment) atoms will either overlap (square inscribed in circle; periodic bound conditions) or atoms will not be included at all in the intensity calculation (cirle inscribed in square). We therefore use square pixels instead of round pixels so that every atom contributes once and only once.

The model must be an integer multiple of the pixel size * sqrt(2) in all directions. There are no error checks for this, you must check this yourself. For models larger than 1 pixel by 1 pixel, the model is segmented into pixel sized chunks and the intensity calculation is run on each pixel, like in the experiment.

There is some common confusion about why the pixel size is what it is. See the pdf included in this repo for more information.

### Alpha
Alpha is the weighting parameter between the fit to experimental data and model energy. This is given on line 12 of the param_file.in file and should be 0.0256, which corresponds to "fit the exerpimental data within the error bars." Ask Jason if you want to know more. This parameter can be fiddled with some. Higher values of alpha indicate fitting the experimental data more heavily.

### Acceptance rate
We ideally want an acceptance rate around 50% throughout the HRMC simulation. This is especially difficult to acheive at low temperatures, and is usually inherently impossible with MC methods unless internally modifying the atom move distance and temperature based upon the acceptance rate. We manually decrease the maximum move distance and temperature throughout the simulation based on the number of steps since the last temperature/max move drop. This number of steps is specified in the parameter file on line 5, after the starting step #.

### Temperature
Setting the initial temperature is important. Too high of a temperature keeps the model in a "liquid" state for too long, wasting computation time. Too low of a starting temperature will not allow the model configuration to escape its local minimum in the potential energy landscape. To estimate a starting temperature, run the HRMC code for a few thousand steps then parse the output to include only the values of "Del-chi" at each step. Put these values into a file, and pass it to calc_starting_temp.py. This will give you a decent estimate of what you want the starting temperature to be.

### Energy
The final energy of the model after an HRMC simulation, which may take 3-6 million steps (attempted atom moves), is up for debate. The important comparison is to the energy of a model that was minimized against only the energy (e.g. via LAMMPS or an energy-only refined MC scheme). The goal is for the energy to be as close to that lowest energy as possible. However, this is unrealistic due to the competition of the energy and fit to exerimental data. We hope the energy will converge within 0.0256eV/atom, which is what alpha is based on. JWHs models were about 0.018eV/atom higher than the perfect (I think). The energy calcuation is not identical to LAMMPS.

### Simulation checks
It is important to check on things as the simulation progesses. You should look at the V(k) for the most recent model file. The V(k) may be printed along with the model, but if not just run the model through femsim (another repo). You should also keep a close eye on the chi_squared_\<jobid>.txt file.

### Autoslice
Autoslice, the slow but much more accurate intensity calculation, is implemented but not tested so I do not know if it works correctly. Hopefully the default settings that are in there now are fine.

##### Note: ZrCuAl2011.eam.alloy is the correct potential file for ZrCuAl.
