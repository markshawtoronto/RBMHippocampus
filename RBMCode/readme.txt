To run the simulation, run RunMe.m. The code for an RBM is in RBoltzmann.m. 
RunMe.m currently includes some degree of temporal blurring over time points in the grid cell input.

To plot some of the output grid cells, run Tester.m. You may have to run Shortcut.m first to generate output from the trained weights. 

It is normal for this code to take several hours to run, and I often ran it overnight. 


Modified versions of the code are included in AlternativeModifications for historical purposes.

The grid cell firing data used as input is contained in 100GridCellsFiringInputData.mat.
	* savedAbstract is a 360000 ms of rat trajectory x 102 matrix. 
	* The first 100 columns are simulated grid cell firing data
	* The last 2 are x and y coordinates if I remember correctly.



After all of the weights are trained, run Shortcut.m to regenerate fantasies of what it thinks it's encoding. The code for that is in Downwardspass.m, and it's basically a shortened modified version of the other RBM file. It should be a little different each time because the RBM bases things on probabilities of firing, but most of the time it should be close. Also it sometimes encodes things you don't really want it to (like using both hidden units to learn the same thing)



Originally by Mark Shaw. 
