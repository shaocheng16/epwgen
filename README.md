# epwgen

This script was made in an attempt to streamline the EPW workflow. It will generate all necessary input files to calculate the electron phonon coupling of a system using the EPW package of Quantum Espresso. It features automatic wannierization window determination and automatic splitting of the phonon calculation into calculations for every irreducible representation of each irreducible q-point.

The script will generate a parent directory called $base_dir which itself contains the following directories:
 - ELB (directory for electron band calculations)
 - PHB (phonon bands)
 - EPM (electron phonon matrix, phonon linewidhts and Eliashberg spectral function)
 - ISO (solving the isotropic Eliashberg equations)

All you need to do is to specify the calculation parameters in the main script and then run the submission scripts called
ep_bands.sh and epw.sh on the cluster in that order.
ep_bands.sh will prepare and submit the calculations for the electron and phonon bands. Because image parallelization does not work for phonon dvscf calculations, this script splits the phonon calculation into multiple calculations; one calculation for each irreducible q-point and additionally one calculation for every irreducible representation of each irreducible q-point. This splitting can also be disabled.

Once the electron and phonon band calculations are finished, check the results and then execute epw.sh.
This will update the EPW input files with the data from the electron and phonon bands calculations and then submit them.

Check the submission scripts for further options/instructions.

Any calculation parameters that are missing or have been "hardcoded" should be easy enough to introduce/change in the input lists called *_in. Everything should still work if you wish to do so.

The calculations need to be run on a cluster with IBM® Platform™ LSF® implementation however modifying it to run on different workload schedulers like SLURM should not be difficult (positions of necessary changes are denoted by "LSF" in the main script). 

For spin orbit coupling you will most likely need an updated postprocessing script ("pp.py").

Solving the Eliashberg equations does not work at the moment.

The script contains a crude but quick and working example. You only need to specify the absolute path to the pseudopotentials on your cluster.
