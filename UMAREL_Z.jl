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


  main="/Users/francovazza/Dropbox/Julia_prog/UHECRpropa/UMAREL/UMAREL/" #..main folder containing UMAREL functions
  #...modules 
  include(string(main,"/constants.jl"))
  include(string(main,"CRadvect_assignZ.jl"))   #...external module with all relevant functions used for the transport of CRs
      
     
      #....input folder and file for the ENZO simulation 
      root0=string("/Users/francovazza/Desktop/data/DATA/RADGAL/LEONARDO/42.5Mpc/snap/")
      root0=string("/Users/francovazza/Desktop/data/DATA/MAKITRA/var9/")
      root_out=string(root0,"/out/")
      snap=[7,7,10,10, 15,15,15,  21,21,21,21]   #...list of ENZO snapshots to be used (they can repeat) 
      zeds=[2.0,1.75,1.5,1.2,1.0,0.7,0.5,0.2,0.1,0.05,0.0]   #...list of redshits associated to snapshots. At each change in redshift, new UHECRs are injected 
      nz=size(zeds)

          
        const dx=292.0 #kpc  resolution of each cell 
        xc=1/(dx)
        scale=dx*1e-3*cmtoMpc
         #.....conversion files from ENZO code units into physical units 
        zfin=0.80   #...final redshift 
        zin=2.0     #...initial redshift
        timeU=Gyr*ustrip(lookback_time(cosmo,1000)::Number)   #...age of the Universe 
        time_tot=Gyr*ustrip(lookback_time(cosmo,zin)::Number)   #.....maximum propagation time (in s), computed based on cosmology and the zin given above 
        Δt=timeU-time_tot  #..initial epoch 
        cdd=2.82e-30#    #...from code density to g/cm^3, physical units; multiply for (1+z)^3 for physical units at a given redshift 
        cv=3848793462.5807     #...from code velocity to cm/s, physical units
        cb=sqrt(cdd*4*pi)*cv    #..from code B-field to G, physical units ; ...multiply for (1+z)^0.5 for physical units at a given redshift 
        
       

    #.....main parameters for this UMAREL run of proagation of UHECRs 
    dsource=0.99 #..only used if UHECRs are selected based on overdensity (in "assign_CR_dens_z"), it is the fraction wrt to the maximum density from which the injection starts 
    mass_source=1e12  #.....only used if UHECRs are selected based on halo catalogs ("assign_CR_halo_z"), it is the minimum total mass of halos considered as source of UHECRs

    #...three possibilities to initialise energy
    #...E_initial=[5e20]  -> single value: all particles are initialised with this energy
    #...E_initial=[18.0,18.5,19.0,19.5...] -> N values: particles are randomly assigned one of these values
    #...E_initial=[-1,17,4]  -> negative value: particles are randomly injected from a E^(-1) distribution, starting from log(17) and for 4 decades in energy
  
    E_initial=[-1,17,4] #  [-1,18,4]# [18.0,18.5,19.0,19.5,20.0,20.5,21.0,21.5]  #....initial energy (in eV) of all injected UHECR. Each entry of E_initial represents an energy bin which can be randomly associated with any injected UHECRS
    Z=1            #....nuclear charge Z=1 proton, Z=2 helium,  Z=7 nitrogen Z=26 iron   Only these are supported (only for them we have loss curves)  
    courant=3.0     #...courant condition for time stepping 
    skip_path=20  #...we write the final path only every a number skip_path steps (to save memory)
    #...boundary of the extracted region within the input simulation
    n=512  #...this is the 1D size of the box which is going to be extracted (1024 is the max possible one)
    i1=1
    i2=i1+n-1
    j1=1
    j2=j1+n-1
    l1=1
    l2=l1+n-1
   
    energy,dEdt=lossesC(Z,main)   #....loading the appropriate tabulated loss function 

  #.....propagation of the particle from t=0 to t corresponding to max_it
  dt=courant*scale/vc
  max_it=convert(Int64,trunc(time_tot/dt))      #...maximum number of iterations
  
  println("going to evolve UHECR for ",max_it," iterations, to cover ",time_tot/Gyr, "Gyr of evolution")
  
  #....this array will store the trajectories and energy evolution of all particles 
    np = 10000 #...total number of UHECR 
    p=Array{Float64}(undef,10,np)
    p.=0.0
    ngen=nz[1]
    npath=convert(Int64,trunc(max_it/skip_path))
    path=Array{Float64}(undef,np,10,npath)
    path.=0.0
    times=Array{Float64}(undef,max_it)  #.....it was npath 
    zed=Array{Float64}(undef,max_it)
  

    #...we use the approximated inversion between time and scale factor as in https://academic.oup.com/mnras/article/505/2/2764/6289943, eq.15-17 
    dz=(zin-zfin)/max_it
#    time_tilde=0.666/(ch*100*sqrt(1-cOmegaM))/Gyr
    time_tilde=(6.519/ch)*(1/sqrt(1-cOmegaM))
    @inbounds for i in 1:max_it #...we define an array of times and redshifts
    zed[i]=zin-i*dz+1e-5       
    times[i]=dt*i/Gyr
    a=(sqrt(cOmegaM/(1-cOmegaM))*sinh((Δt/Gyr+times[i])/time_tilde))^0.6666         
    zed[i]=1/a-1     
    end 
  
     iz0=0
     inj=0
     inj0=0
     @inbounds for t in eachindex(times)  #...main time loop 
      global iz0
      global inj, np,inj0
      global path,p

      println("iteration ",t, "of ",max_it, "z=",zed[t])
  
      iz=findz(zeds,zed[t])  #...from the list of available snapshots, we find the one closest to the epoch of this t-iteration
  
         if iz[1] != iz0   #...anytime we switch to the next redshift snapshot, we inject new UHECRs 
         global  inj+=1
          println("NEW INJECTION at z=",zed[t])
         end 
         iz0=iz[1]
  
#    file1=string(root0,"minus1.0_0.2_j60fe3I__full_dtb_R00",snap[iz])
#    file1=string(root0,"minus1.0_0.2_j60fe3I__full_dtb_R0027")
    
      file1=string(root0,"minus2.9_08RD0_dtb_",snap[iz])   #...ENZO snapshots with 3D data 
      filecat=string(root0,snap[iz],"_minus2.9_08_halof_100_new.dat")   #...halo catalogs 
      if snap[iz] <10
      file1=string(root0,"minus2.9_08RD0_dtb_0",snap[iz])
      filecat=string(root0,"0",snap[iz],"_minus2.9_08_halof_100_new.dat")
      end 
      if inj>inj0
      inj0=inj
     # p=assign_CR_dens_z(p,i1,i2,j1,j2,l1,l2,file1,dsource,E_initial,np,inj,ngen)   #...assign UHECR based on density 
      p=assign_CR_halo_z(p,i1,i2,j1,j2,l1,l2,filecat,mass_source,E_initial,np,inj,ngen,n)   #...assign UHECR based on density  
      inj0=inj
      ntr=size(p)    #....number of UHECR - it depends on the dsource density choosen above 
      np=ntr[2]
      npt=convert(Int64,trunc(np/ngen))   #....number of new UHECR to be injected at this step 
      nj1=(inj-1)*npt+1
      nj2=nj1+npt-1
     #....we plot the initial location of UHECR
      plo=plot(p[1,nj1:nj2],p[2,nj1:nj2],seriestype=:scatter,ms=0.01,label="",grid=false,aspect_ratio=1.0)  #...plot initial location of CRs
      xlims!(i1, i2)
      ylims!(j1,j2)
      #zlims!(l1,l2)   #...if this is uncommented, it is a 3D plot 
      filep1=string(root_out,"_initial_map_newU",snap[iz],".png")
      savefig(filep1)

    bx=h5read(file1,"Bx",(i1:i2,j1:j2,l1:l2))
    by=h5read(file1,"By",(i1:i2,j1:j2,l1:l2))
    bz=h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2))

  global    bx*=cb*(1+zed[t])^2    #...this field is comoving. it is turned into physical, x(1+z)^2, inside the propagation routine
  global    by*=cb*(1+zed[t])^2
  global    bz*=cb*(1+zed[t])^2
  end 
#....main function which does all the hard work 

  path=move_CR(p,t,skip_path,scale/(1+zed[t]),courant,dt,dx,i1,i2,j1,j2,l1,l2,path,cdd,cv,cb,energy,dEdt,Z,zed[t],ngen,inj,file1,bx,by,bz) #...evolve CR in time 

     end 

  #....plotting and writing of data on disk 
  
          #1) static plot without boundary conditions
            xs=dx*1e-3  #...to convert positions from cells to Mpc
            jump=20
            #the following plot can use a lot of memory, do for i in 1:jump:np with jump = any integer number to plot only 1 every jump particles
            @inbounds for i in 1:jump:np-1   #....all particle trajectories are plotted without periodic BC 
            if i==1 
            plo=plot(path[i,5,:]*xs,path[i,6,:]*xs,label="",aspect_ratio = :equal ,dpi=1000,lw=0.0,ms=0.3,grid=false,seriestype=:scatter)
            end 
            plot!(path[i,5,:]*xs,path[i,6,:]*xs,label="",lw=0.0,alpha=0.5,seriestype=:scatter,ms=0.1)
          end
            tit=string(" E0=",E_initial)

            if E_initial[1]==-1
              tit=L" N(E) \propto E^{-1}"
            end 

            title!(string("Z=",Z," ,",tit," eV"),fonts=20)
                yaxis!("[Mpc]",fonts=20)
           xaxis!("[Mpc]",fonts=20)
           filep1=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_noBC.png")
           savefig(filep1)

          #2) animated gif with periodic boundary conditions
        #...plotting evolution of UHECR in the periodic BC as an animated gif 
        it=0
                anim=@animate for t in 1:skip_path:npath*skip_path   #....plotting only the first 100 steps (soon after that all is isotropic! )
            global      it+=1
        println(t)      
        tit=string(" E0=",E_initial)
        if E_initial[1]==-1
          tit=L" N(E) \propto E^{-1}"
        end 
        npt=convert(Int64,trunc(np/ngen)) 
        plot(path[1:2,1,it]*xs,path[1:2,2,it]*xs,label="",aspect_ratio = :equal,line=0,dpi=300,lw=0.0,marker=:circle,ms=0.1,grid=false)
 #  
    cg=["red","blue","green","brown","cyan","blue","pink","grey","yellow","orange","magenta","green","blue","red","pink","orange","gold"]
@inbounds for gg in 1:nz[1]
          g1=(gg-1)*npt+1
          g2=g1+npt-1
  #        g1=1
   #       g2=

        plot!(path[g1:g2,1,it]*xs,path[g1:g2,2,it]*xs,label="",line=0,lw=0.0,marker=:circle,ms=1.0,seriestype=:scatter,color=cg[gg],alpha=0.5)
        end 
        title!(string("Z=",Z," ,",tit," eV, time= ",trunc(times[it],digits=4),"Gyr"),fonts=20)
        yaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        xaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        end 
        file_gif=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_BC_redshift_newU.gif")
        gif(anim,file_gif,fps=10)    
        


        #...for safety, HDF5 files must be manually deleted before overwriting - otherwise the code stops here
        filep1=string(root_out,"path_",E_initial[1],"Z",Z,"_spec_redshift.hdf5")
        h5write(filep1,"px_pbc",path[:,1,:])  #...X position with periodic BC 
        h5write(filep1,"py_pbc",path[:,2,:])  #...Y position with periodic BC 
        h5write(filep1,"pz_pbc",path[:,3,:])  #...Z position with periodic BC 
        h5write(filep1,"px",path[:,5,:])  #...X position without periodic BC 
        h5write(filep1,"py",path[:,6,:])  #...Y position without periodic BC 
        h5write(filep1,"pz",path[:,7,:])  #...Z position without periodic BC 
        h5write(filep1,"E[eV]",path[:,4,:])  #...energy 
        h5write(filep1,"B[G]",path[:,8,:])  #...physical mag.field in Gauss
        h5write(filep1,"redshift",path[:,9,:])  #...redshift
        h5write(filep1,"time",path[:,10,:])  #...elapsed time  

#        h5write(filep1,"time",times)         #...list of simulated timesteps
    ##   h5write(filep1,"B[G]",path[:,8,:])  #...comoving mag.field in Gauss
       
        println("end of run")
