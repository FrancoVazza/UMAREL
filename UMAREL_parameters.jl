#...parameters of the simulation

  #...COSMOLOGY
  using Cosmology
  global cOmegaM=0.308
  global ch=0.678
  global cosmo=cosmology(OmegaM=cOmegaM,h=ch)

      #....input folder and file for the ENZO simulation 
      root0=string("/Users/francovazza/Desktop/data/DATA/MAKITRA/var11/snap/")
      root_halos=string("/Users/francovazza/Dropbox/DATA/MAKITRA/var11/")
      root_out=string(root0,"/out/")
      snap=[13, 14, 15, 17, 19, 19,  21, 21,21]   #...list of ENZO snapshots to be used (they can repeat) 
      zeds=[1.84,1.455,1.065,0.555,0.3,0.202,0.05,0.02,0.0]   #...list of redshits associated to snapshots. At each change in redshift, new UHECRs are injected 
      #snap=[19, 21,21]   #...list of ENZO snapshots to be used (they can repeat) 
      #zeds=[0.05,0.02,0.0]   #...list of redshits associated to snapshots. At each change in redshift, new UHECRs are injected 
    
      nz=size(zeds)
     

       #..VARIOUS PARAMETERS
     
        const dx=292.0 #kpc  comoving resolution of each cell 
        xc=1/(dx)     #...useful 
        scale=dx*1e-3*cmtoMpc
         #.....conversion files from ENZO code units into physical units 
        zfin=0.0   #...final redshift 
        zin=zeds[1]
        timeU=Gyr*ustrip(lookback_time(cosmo,1000)::Number)   #...age of the Universe 
        time_tot=Gyr*ustrip(lookback_time(cosmo,zin)::Number)-Gyr*ustrip(lookback_time(cosmo,zfin)::Number)   #.....maximum propagation time (in s), computed based on cosmology and the zin given above 
        Δt=timeU-time_tot  #..initial epoch 
        cdd=2.82e-30#    #...from code density to g/cm^3, physical units; multiply for (1+z)^3 for physical units at a given redshift 
        cv=3848793462.5807     #...from code velocity to cm/s, physical units
        cb=sqrt(cdd*4*pi)*cv    #..from code B-field to G, physical units ; ...multiply for (1+z)^0.5 for physical units at a given redshift 
        
       

      #.....main parameters for this UMAREL run of proagation of UHECRs 
      dsource=0.99 #..only used if UHECRs are selected based on overdensity (in "assign_CR_dens_z"), it is the fraction wrt to the maximum density from which the injection starts 
      mass_source=1e13  #.....only used if UHECRs are selected based on halo catalogs ("assign_CR_halo_z"), it is the minimum total mass of halos considered as source of UHECRs
      #...three possibilities to initialise energy
      #...E_initial=[5e20]  -> single value: all particles are initialised with this energy
      #...E_initial=[18.0,18.5,19.0,19.5...] -> N values: particles are randomly assigned one of these values
      #...E_initial=[-1,17,4]  -> negative value: particles are randomly injected from a E^(-1) distribution, starting from log(17) and for 4 decades in energy
      E_initial=[-1,17,4] #  [-1,18,4]# [18.0,18.5,19.0,19.5,20.0,20.5,21.0,21.5]  #....initial energy (in eV) of all injected UHECR. Each entry of E_initial represents an energy bin which can be randomly associated with any injected UHECRS
      Z=1            #....nuclear charge Z=1 proton, Z=2 helium,  Z=7 nitrogen Z=26 iron   Only these are supported (only for them we have loss curves)  
      courant=90.0     #...courant condition for time stepping 
      skip_path=20  #...we write the final path only every a number skip_path steps (to save memory)
      cosmo=1     #...if cosmo=1 use cosmology expansion factors, 1+z dependences etc ;  if cosmo=0 we still use the list of redshift given above, but without z-dependent effect 
      use_syntheticB=0  #....0=using cosmological simulations, 1=using a test 3D field 
      #...boundary of the extracted region within the input simulation
      n=512  #...this is the 1D size of the box which is going to be extracted 
     global  i1=1
     global i2=i1+n-1
     global j1=1
     global j2=j1+n-1
     global l1=1
     global l2=l1+n-1
   
    