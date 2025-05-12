

  using LaTeXStrings
  using Plots
  Plots.PyPlotBackend()
  using Statistics
  using SpecialFunctions
  using DelimitedFiles
  using HDF5
  using Unitful
  using Random

#...cosmological parameters
using Cosmology
cOmegaM=0.308
ch=0.678
cosmo=cosmology(OmegaM=cOmegaM,h=ch)


  main="/Users/francovazza/Dropbox/Julia_prog/UHECRpropa/MAGHOR/MAGHOR/" #..main folder containing ROGER functions
  #...modules 
  include(string(main,"/constants.jl"))
  include(string(main,"CRadvect_assign.jl"))   #...external module with all relevant functions used for the transport of CRs
      

    #......selection of input cosmology and files 
    cosmo=cosmology(OmegaM=0.258,h=0.72)   #...setting cosmology here
     
      #....input folder and file for the ENZO simulation 
      root0=string("/Users/francovazza/Desktop/data/DATA/RADGAL/LEONARDO/42.5Mpc/snap/")
      root_out=string(root0,"/out/")
      file1=string(root0,"minus1.0_0.2_j60fe3I__full_dtb_R0027")
      const dx=41.5 #kpc  resolution of each cell 
      xc=1/(dx)
      scale=dx*1e-3*cmtoMpc
       #.....conversion files from ENZO code units into physical units 
       z=0.02   #...redshift 
       cdd=2.82e-30*(1+z)^3. #...from code density to g/cm^3, physical units; divide for 1/(1+z)^3 for comoving
       cv=1.258e+08 #...from code velocity to cm/s, physical units; multiply this for 1+z for comoving
       cb=sqrt(cdd*4*pi)*cv# #..from code B-field to G, physical units ; ...divide for 1/(1+z)^2 for comoving
       t0=ustrip(lookback_time(cosmo,z)::Number)
      #...............................................
      


    #.....main parameters for this MAGHOR run of proagation of UHECRs 
    dsource=250 #...lower gas density threshold for the injection of UHECRs, relative to the cosmic mean gas density. 
    E_initial=[-1]# [18.0,18.5,19.0,19.5,20.0,20.5,21.0,21.5]  #....initial energy (in eV) of all injected UHECR. Each entry of E_initial represents an energy bin which can be randomly associated with any injected UHECRS
    Z=1            #....nuclear charge Z=1 proton, Z=2 helium,  Z=7 nitrogen Z=26 iron   Only these are supported (only for them we have loss curves) 
    time_tot=9e16   #.....maximum propagation time (in s)  (3e16->1Gyr)
    courant=9.0     #...courant condition for time stepping 
    skip_path=5    #...we write the final path only every a number skip_path steps (to save memory)
    #...boundary of the extracted region within the input simulation
    n=300   #...this is the 1D size of the box which is going to be extracted (1024 is the max possible one)
    i1=100
    i2=i1+n-1
    j1=100
    j2=j1+n-1
    l1=200
    l2=l1+n-1
   
    energy,dEdt=losses(Z,main,z)   #....loading the appropriate tabulated loss function 

  println("Injecting new UHECR based on density")
  trac=assign_CR_dens(i1,i2,j1,j2,l1,l2,file1,dsource,E_initial)   #...assign UHECR based on density 
  ntr=size(trac)    #....number of UHECR - it depends on the dsource density choosen above 
  np=ntr[2]
  println("A total number of ",np, " UHECR is going to be simulated")


  #....we plot the initial location of UHECR
  plo=plot(trac[1,:],trac[2,:],seriestype=:scatter,ms=0.01,label="",grid=false,aspect_ratio=1.0)  #...plot initial location of CRs
  xlims!(i1, i2)
  ylims!(j1,j2)
  #zlims!(l1,l2)   #...if this is uncommented, it is a 3D plot 
  filep1=string(root_out,"_initial_map.png")
  savefig(filep1)
  
  #.....propagation of the particle from t=0 to t corresponding to max_it
  dt=courant*scale/vc
  max_it=convert(Int64,trunc(time_tot/dt))      #...maximum number of iterations
  println("going to evolve UHECR for ",max_it," iterations, to cover ",time_tot/3e16, "Gyr of evolution")
  
  #....this array will store the trajectories and energy evolution of all particles 
  npath=convert(Int64,trunc(max_it/skip_path))
  path=Array{Float64}(undef,np,7,npath)
   path.=0.0
   times=Array{Float64}(undef,npath)
   for t in eachindex(times)  #...convenient array of times (in Gyr)
    times[t]=dt*t*skip_path/(1e9*yr)
    end 
   #....main function which does all the hard work 
  path=move_CR(trac,max_it,skip_path,scale,courant,dt,dx,i1,i2,j1,j2,l1,l2,path,cdd,cv,cb,energy,dEdt,Z) #...evolve CR in time 

   println("evolution done")

  #....plotting and writing of data on disk 
  
          #1) static plot without boundary conditions
            xs=dx*1e-3  #...to convert positions from cells to Mpc
            jump=200
            #the following plot can use a lot of memory, do for i in 1:jump:np with jump = any integer number to plot only 1 every jump particles
            @inbounds for i in 1:jump:np   #....all particle trajectories are plotted without periodic BC 
            if i==1 
            plo=plot(path[i,5,:]*xs,path[i,6,:]*xs,label="",aspect_ratio = :equal ,dpi=1000,lw=0.1,grid=false)
            end 
            plot!(path[i,5,:]*xs,path[i,6,:]*xs,label="",lw=0.2,alpha=0.2)
            end 
           title!(string("Z=",Z," E0=",E_initial," eV"),fonts=20)
           yaxis!("[Mpc]",fonts=20)
           xaxis!("[Mpc]",fonts=20)
           filep1=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_noBC.png")
           savefig(filep1)
    
          #2) animated gif with periodic boundary conditions
        #...plotting evolution of UHECR in the periodic BC as an animated gif 
        nmaps=50
        step=convert(Int64,trunc(max_it/nmaps))
        println(max_it," ",step)
         # anim=@animate for t in 1:step:max_it   #....plotting only every "step" timestep, to produce a total of 100 maps 
          anim=@animate for t in 1:100   #....plotting only the first 100 steps (soon after that all is isotropic! )
          
        println(t)      
        plot(path[:,1,t]*xs,path[:,2,t]*xs,label="",aspect_ratio = :equal ,line=0,dpi=300,lw=0.0,marker=:circle,markercolor=:blue,color="white",ms=0.1,grid=false)
        title!(string("Z=",Z," E0=",E_initial," eV, time= ",trunc(times[t],digits=4),"Gyr"),fonts=20)
        yaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        xaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        end 
        file_gif=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_BC.gif")
        gif(anim,file_gif,fps=10)    
        


        #...for safety, HDF5 files must be manually deleted before overwriting - otherwise the code stops here
        filep1=string(root_out,"path_",E_initial[1],"Z",Z,"_spec.hdf5")
        #h5write(filep1,"px",path[:,1,:])  #...X position with periodic BC 
        #h5write(filep1,"py",path[:,2,:])  #...Y position with periodic BC 
        #h5write(filep1,"pz",path[:,3,:])  #...Z position with periodic BC 
        h5write(filep1,"px",path[:,5,:])  #...X position without periodic BC 
        h5write(filep1,"py",path[:,6,:])  #...Y position without periodic BC 
        h5write(filep1,"pz",path[:,7,:])  #...Z position without periodic BC 
        h5write(filep1,"E[eV]",path[:,4,:])  #...energy 
        h5write(filep1,"time",times)         #...list of simulated timesteps
       
        println("end of run")