  using LaTeXStrings
  using Plots
  Plots.PyPlotBackend()
  using Statistics
  using SpecialFunctions
  using DelimitedFiles
  using HDF5
  using Unitful
  using Random

  #...just two parameters to give here 
  tag="B0_run"    #....will be attached to all output file names to differentiate them if necessary 
  np = 20000   #...total number of UHECR     
  
    #..PATH TO FOLDERS OF ROUTINES AND FILES
      main="/Users/francovazza/Dropbox/Julia_prog/UHECRpropa/UMAREL/UMAREL/" #..main folder containing UMAREL functions
      include(string(main,"/constants.jl"))
      include(string(main,"parameters_UMAREL.jl"))
      include(string(main,"functions_UMAREL.jl"))   #...external module with all relevant functions used for the transport of CRs    
      energy,dEdt=lossesC(Z,main)   #....loading the appropriate tabulated loss function 
      times,zed,dt,max_it=define_times(courant,scale,zfin,zin,ch,cOmegaM)   #...computes the maximum iterations and the array of redshift to cover the entire evolution
   
    
     #....this array will store the trajectories and energy evolution of all particles 
     p=Array{Float64}(undef,10,np)
     p.=0.0
     ngen=nz[1]
     npath=convert(Int64,trunc(max_it/skip_path))
     path=Array{Float64}(undef,np,10,npath)
     path.=0.0
    

    #....MAIN SIMULATION CODE 

     iz0=0   #...useful counters 
     inj=0
     inj0=0
     @inbounds for t in eachindex(times)  #...main time loop 
     global iz0, inj, np,inj0, path,p,cosmo,use_syntheticB

      println("iteration ",t, "of ",max_it, "z=",zed[t])
   
       iz=findz(zeds,zed[t])  #...from the list of available snapshots, we find the one closest to the epoch of this t-iteration
    
       #...INJECTION OF NEW UHECR 
        if iz[1] != iz0   #...anytime we switch to the next redshift snapshot, we inject new UHECRs 
            global inj+=1
            if inj>inj0   
            inj0=inj    
            println(root0,root_halos,snap[iz],i1,i2,j1,j2,l1,l2,mass_source,E_initial,np," ",inj," ",ngen," ",n)
            global p,bx,by,bz= inject_new_UHECR(root0,root_halos,snap[iz],p,i1,i2,j1,j2,l1,l2,mass_source,E_initial,np,inj,ngen,n)
            end     
            println("NEW INJECTION at z=",zed[t])
        end 
         iz0=iz[1]

  

      #...MAIN PROPAGATION & LOSSES ROUTINE
     path=move_CR(p,t,skip_path,scale,courant,dt,dx,i1,i2,j1,j2,l1,l2,path,cdd,cv,cb,energy,dEdt,Z,zed[t],ngen,inj,bx,by,bz) #...evolve CR in time 

     end 

   println("the propagation of UHECRs is done, now plotting")
    
  #....plotting and writing of data on disk 
  
       a=plot_map_WBC(dx,np,path,E_initial,Z,root_out,cosmo,tag)   #...plot trajectories without periodic BC 
       a=make_gif(skip_path,npath,E_initial,path,nz,Z,root_out,n,cosmo,tag)    #...gif with UHECR propagation
       a=plot_spectrum(path,np,root_out,E_initial,Z,cosmo,tag)
       a=plot_spectrum_suppression(path,np,root_out,E_initial,Z,cosmo,tag,root_halos)
       a=write_path(root_out,E_initial,Z,cosmo,tag,path)


        println("end of run")
