# MAGHOR
 code for the propagation of UHECRs on ENZO simulations
 (work in progress)


These are free parameters which the user should set to start a run:

    
* dsource=200     #...lower gas density threshold for the injection of UHECRs, relative to the cosmic mean gas density. 
* E_initial=1e18   #....initial energy (in eV) of all injected UHECR
* Z=1              #....nuclear charge Z=1 proton, Z=2 helium,  Z=7 nitrogen Z=26 iron   Only these are supported (only for them we have loss curves) 
* time_tot=3e16    #.....maximum propagation time (in s)  (3e16->1Gyr)
* courant=3.0      #...courant condition for time stepping 
* n=400   #...this is the 1D size of the box which is going to be extracted (1024 is the max possible one)


 This is a movie showing the propagation of injected cosmic rays for the first 100 timesteps
<img src="_UHECR_path_color_map_1.0e18_Z_26_BC.gif" alt="alt text" width="whatever" height="whatever">


This is the full set of trajectories, unfolding the effect of periodic boundary conditions
<img src="_UHECR_path_color_map_1.0e18_Z_26_noBC.png" alt="alt text" width="whatever" height="whatever">
