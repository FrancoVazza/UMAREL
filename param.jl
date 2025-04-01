#...parameters of the simulation
@everywhere const n=256 #grid size
@everywhere const dx=18.04   #kpc
@everywhere xc=1/(dx)
@everywhere run_model=["MH1"]                                                     #....list of runs to analyse (just one for the moment )
@everywhere rrr=1
@everywhere runm=run_model[rrr]


@everywhere rootf=string("/Users/francovazza/Desktop/data/DATA/MH/",runm,"/")     #..folder containing the hdf5 files of tracers
@everywhere root_map=string("/Users/francovazza/Desktop/data/DATA/MH/",runm,"/")  #..folder where the 2D maps will be put
@everywhere root_out=string("/Users/francovazza/Desktop/data/DATA/MH/",runm,"/")  #..foldere where the spectra will be put
@everywhere snap_in0=[7]  #..initial snapshot
@everywhere ntr00=[106442] #...maximum possible number of tracers
@everywhere snap_fin0=[115] #...final snapshot
@everywhere ntr0=ntr00[rrr]
@everywhere snap_in=snap_in0[rrr]
@everywhere snap_fin=snap_fin0[rrr]
@everywhere initial_snap=snap_in

#...cosmological parameters
@everywhere using Cosmology
@everywhere cOmegaM=0.308
@everywhere ch=0.678
@everywhere cosmo=cosmology(OmegaM=cOmegaM,h=ch)


#...parameters for spectra
@everywhere using SpecialFunctions
@everywhere const  p_max=log10(5e6)
@everywhere const  p_min=log10(1.0)
@everywhere const  dp=0.1  # log spacing of momenta in the spectra
@everywhere const  part=2  #...1=proton  2=electron
@everywhere const np=1+(floor(Int64,(p_max-p_min)/dp))
@everywhere const pend=np
@everywhere pval=collect(p_min:dp:p_max)
# parameters for shock modelling
@everywhere const norm=30 #to have all in 1e30 erg or 1e30 erg/s
@everywhere const lcorr=20.0 #...in kpc
@everywhere const mthr=1.7   #...minimum Mach number for shock re-acceleration (i.e. no injection of fresh new electrons)
@everywhere const minj=2.5   #...minimum Mach number for the shock injection (new fresh electrons are added to spectra)
#..other parameters
@everywhere const ηB=0.05   #...dissipation efficiency of turbulence into B-fields
@everywhere Dmin=5e-27   #...minimum physical density to start the injection of tracers in the run
@everywhere const scale_turbo=2*dx   #....stencil used to compute vorticity


#....radio frequencies used to compute synchrotron emission
@everywhere freqf=[5e7,1.5e8,6.5e8,1.4e9,5e9]#
#@everywhere freqf=[1e6,3e6,1e7,3e7,1e8,3e8,1e9,3e9,1e10]#
@everywhere nfreq=size(freqf)
