*Modifications*

These scripts were modified by Mark Shaw in 2012 to output the abstract grid cell firing data only and heavily commented. The rat trajectory comes from Hafting et al. (2005).

***********

Grid cell network simulations. Version 1.02.

These MATLAB (R14) scripts are cleaned up versions of the scripts used
to generate the manuscript figures. These will allow interested
readers to verify the results or to easily continue this line of
research by examining further directions. Note: Most Figure numbers
were changed during revisions--while I have tried to update the figure
references in these scripts, it is possible I have missed one or more.
Eric Zilli 2010 June 1

The scripts use print_eps to save the figures as an eps file. This is
available on MATLAB Central: It is really useful for making figures!
(v1.02 note: print_eps is now called export_fig but is just as good)

This archive includes experimental rat trajectory information in the
file rat_10925.mat. This trajectory freely available on the Moser lab
website and was collected for the Hafting et al 2005 paper by members
of the Moser lab. Cheers to them for making this information freely
available.

The general process resulting in a 2D grid cell simulation given some
model of an individual cell is:

1. Find injected current levels (variable I in the paper, I and
inputMags in the scripts) and noise level (sigma in the paper,
uniqueNoiseSTD in the scripts) that match the biological median ISI of
0.2 s period (5 Hz firing rate) and a standard deviation of 0.036
s. Because noise will cause each run to have a different mean and
std. dev. of the ISIs, it is useful to simulate, say, 5000 uncoupled
cells (allncells=5000, pcon=0) and check the statistics of the median
cell.

This can be done with scripts:
* SI_simple_model_stability_vs_params.m
* SI_acker_model_stability_vs_params.m

The output of the script is a string of numbers (parameter values and
measures of the network's activity) which are optionally written to a
file (set fileOut=1). Search for fprintf for the meaning of the values.

After running a simulation, try (some variables may not be produced by
all scripts):
% voltage variable and gating/auxilary variables
% for cell #1 during entire simulation (don't forget to transpose!):
figure; plot(state'); 
% voltage variable only:
figure; plot(state(1,:));
% activity of postsynaptic cell during entire simulation
figure; plot(post); 
% mean ISI for each individual cell in the network
figure; plot(indivmeans); 
% ISI standard deviations for each individual cell in the network
figure; plot(indivstds); 
% estimated grid stability time for each individual cell in the
  network
figure; plot(indivstabilities); 

2. Find input-frequency relation (FI curve). This should be done to a
fairly high resolution (but low resolutions haven't been thoroughly
explored, so perhaps there is a coarser resolution that still works
fine), so the following guidelines may be useful in maximizing the use
of simulation time by simulating no more than is necessary: Use the
desired spacing of the grid cell to find the beta parameter (beta =
sqrt(3)*spacing/2). Now select a maximum instantaneous velocity to be
able to accurately path integrate (say v_max = 1 m/s). The range of
frequencies the FI curve then needs is 2*beta*v_max wide. Next select
a low frequency bound and find the input current that produces it
(I_low).  This is essentially arbitrary for a Class 1 excitable cell,
but for a Class 2 excitable cell the cell itself will impose a minimum
firing rate (by definition). Find the increment (dI) to the input
current that moves it, say 1/200th (a reasonably high resolution) of
the way to the desired upper frequency (which is freq_high = freq_low
+ 2*beta*v_max corresponding to I_high). Now run FI simulations for
inputMags=I_low:dI:I_high.

This can be done with scripts:
* SI_simple_model_FI_relation.m
* SI_acker_model_FI_relation.m

We include pre-generated FI curves for the cases included in the
manuscript. These are:
Acker_sn_FI_n250.mat                  Synaptically-coupled network of 250
                                      noisy biophysical cells
simple_model_RS1_FI_Jan09_n1.mat      Single simple model, Type 1
                                      excitable, regular spiking cell
simple_model_RS1n_FI_Jan09_n1.mat     Noisy single simple model, Type 1
                                      excitable, regular spiking cell
simple_model_RS2_FI_Jan09_n1.mat      Single simple model, Type 2
                                      excitable, regular spiking cell
simple_model_RS2n_FI_Jan09_n1.mat     Noisy single simple model, Type 2
                                      excitable, regular spiking cell
simple_model_RS1gn_FI_Jan09_n250.mat  Gap-junction--coupled network of
                                      250 noisy, Type 1 regular spiking
                                      simple model cells
simple_model_RS1sn_FI_Jan09_n250.mat  Synaptically-coupled network of 250
                                      noisy, Type 1 regular spiking
                                      simple model cells
simple_model_RS2gn_FI_Jan09_n250.mat  Gap-junction--coupled network of
                                      250 noisy, Type 2 regular spiking
                                      simple model cells
simple_model_RS2sn_FI_Jan09_n250.mat  Synaptically-coupled network of 250
                                      noisy, Type 2 regular spiking
                                      simple model cells

These files contain vectors "currents" and "freqs". currents(i) is the
injected current level needed to drive a network to fire at freqs(i)
Hz.

3. Use the FI curve in a 2D grid simulation to translate the velocity
signals into desired frequencies. This allows the VCOs to be
controlled, but you still need to find appropriate parameters for the
postsynaptic cell (G in the paper). I have not solved this problem. A
reasonable start is to run 4 s simulations and compare traces of VCO
cells to the activity in the postsynaptic cell (like the figures of
traces in the manuscript) and visually decide whether stronger weights
or time constants (or damping constants for a resonant postsynaptic
cell) need to change and by how much. Once apparently successful
parameters are found, try running a 10 s simulation then a 40 s
simulation then a 180 s simulation (to give approximate magnitudes),
making changes to parameters at any stage where it becomes clear you
have them wrong. This can take a while. I have one analytical and one
numerical technique that are starting points for successfully handling
this problem, but neither proved perfect in practice so are not
included.

This can be done with scripts:
* SI_simple_model_2d_grid.m
* SI_acker_model_2d_grid.m

---

The manuscript also examined firing rate adaptation. This aspect of
the problem can be explored with the scripts:
* SI_simple_model_FI_history_dependence.m
* SI_acker_model_FI_history_dependence.m

---

Parameters in some scripts depend on others:

2d_grid and FI_history_dependence scripts with single noiseless cells
depend on:
* FI_relation for the frequency-current relations needed to control
the single-cell VCOs

2d_grid and FI_history_dependence scripts with noisy cells or networks
depend on:
* stability_vs_params for values of uniqueNoiseSTD (to match biology)
  and then ncells, g, and pcon to achieve high stability times despite
  noise
* FI_relation for the frequency-current relations needed to control
  the network VCOs
FI_relation with noisy cells or networks depends on:
* stability_vs_params for values of uniqueNoiseSTD (to match biology)
  and then ncells, g, and pcon to achieve high stability times despite
  noise

MANUSCRIPT ERRATA

 None yet.

CHANGELOG

New in 1.02 (2011Jan10):

 * Fixed plotType reference to manuscript figure numbers in comments
 in SI_simple_model_2d_grid.m

 * Fixed mistake in hafting_trajectory.m which was not filtering the
 trajectory at the same frequency as in the manuscript.

 * hafting_trajectory no longer plots the filtered and unfiltered
 trajectories (uncomment those lines in hafting_trajectory.m if those
 plots are desired)

 * In SI_simple_model_2d_grid.m simulation type comments, replaced
 variable "D" with "noise".

 * In SI_simple_model_2d_grid.m simulation type 2 has much nicer
 looking parameters for the resonate-and-fire model. Type 7 has nicer
 looking parameters for all postsynaptic models.

New in 1.01:
 * Included the trajectory from Hafting et al. 2005 which was used to
 make the manuscript figures.
