  
function cross_prod(a,b)
   c=similar(a)
   c[1] = a[2]*b[3]-a[3]*b[2]
   c[2] = a[3]*b[1]-a[1]*b[3]
   c[3] = a[1]*b[2]-a[2]*b[1]
   return c 
end

function findz(a,x)
na=size(a)
i0=1#na[1]
   @inbounds for i in eachindex(a)
      if i==1 
      continue
      end 
   if x>=a[i] && x<a[i-1] 
     i0=i
   end 
  end 

 return i0
end 

  function assign_CR_dens_z(p,i1,i2,j1,j2,l1,l2,file1,dCR,E_initial,np,inj,ngen)
    #...this function assigns the initial position of UHECRs within cells denser than something ("halos")
    #...depending on the dCR density threshold, more or less UHECRs are injected 
#    d=h5read(file1,"Dark_Matter_Density",(i1:i2,j1:j2,l1:l2)) 
    d=  h5read(file1,"Grid00000001/Density",(i1:i2,j1:j2,l1:l2)) 
    mad=maximum(d)


    i_source=findall(d.>dCR*mad)
    nd=size(i_source)
#    println("# of UHECR sources=", nd[1])
    npt=convert(Int64,trunc(np/ngen))   #....number of new UHECR to be injected at this step 


    @inbounds for i in 1:npt  #...loop over the number of new UHECRs 
      nE=size(E_initial)
      ip=npt*(inj-1)+i        #...counter to begin the generation starting from the last UHECR evolved this far       
      ir=convert(Int64,trunc(nd[1]*rand()+1))  #...selects a random cell from the i_source selection done above
      ijk = CartesianIndices(d)[i_source[ir]]
            ix=ijk[1]+rand()
            iy=ijk[2]+rand()
            iz=ijk[3]+rand()
            θ=pi*rand()
            φ=2*pi*rand()        
            vx=vc*sin(θ)*cos(φ)     #...random initial velocity vector (v=c)
            vy=vc*sin(θ)*sin(φ)
            vz=vc*cos(θ)
         p[1:3,ip].=[ix,iy,iz]
         p[4:6,ip].=[vx,vy,vz]
         p[7:9,ip].=[ix,iy,iz]

         if nE[1]==1 && E_initial[1]>0
         p[10,ip]=10^E_initial[1]  #...initial energy in eV 
         end 
         if nE[1]>1
            ii=convert(Int64,round(nE[1]*rand()))
            if ii<1
            ii=1
            end 
            if ii>nE[1]
            ii=nE[1]
            end 
            p[10,ip]=10^(E_initial[ii])  #....this generates are random distribution of initial energies picked from the given energy bins.
            end 

            if nE[1]==1 && E_initial[1]==-1
               Er=4*rand()
               Ein_random=Er[1]
               p[10,ip]=10^(Ein_random)  #....this generates are random distribution of initial energies picked from the given energy bins.
               end 
         
         end 
    
    d=nothing
   return p
   end
 

  function assign_CR_halo_z(p,i1,i2,j1,j2,l1,l2,filecat,mass_source,E_initial,np,inj,ngen,n)
   #...this function assigns the initial position of UHECRs using the position of matter halos given by a catalog

   a=readdlm(filecat)

   imass=findall(a[:,5].>=mass_source)
     mass=a[imass,5]
     radius=a[imass,4]
     x=a[imass,1]
     y=a[imass,2]
     z=a[imass,3]
     x.*=n
     y.*=n
     z.*=n
     x.+=(-i1+1)
     y.+=(-j1+1)
     z.+=(-l1+1)
     
     im=sortperm(mass,rev=true)    
  
     npt=convert(Int64,trunc(np/ngen))   #....number of new UHECR to be injected at this step 
     nim=size(im)
     println("sources with >Mass= ", nim, " ",filecat)

     if nim[1]<=1
      imass=findall(a[:,5].>=mass_source*0.1)
      mass=a[imass,5]
      radius=a[imass,4]
      x=a[imass,1]
      y=a[imass,2]
      z=a[imass,3]
      x.*=n
      y.*=n
      z.*=n
      x.+=(-i1+1)
      y.+=(-j1+1)
      z.+=(-l1+1)
      im=sortperm(mass,rev=true) 
      nim=size(im)
      println("sources with >Mass= ", nim, " ",filecat)
 
     end 

     @inbounds for i in 1:npt  #...loop over the number of new UHECRs 
     nE=size(E_initial)
     ip=npt*(inj-1)+i        #...counter to begin the generation starting from the last UHECR evolved this far       
   if i<=nim[1]
      x1=x[im[i]]
      y1=y[im[i]]
      z1=z[im[i]]
   end 
      if i>nim[1]
       ni=convert(Int64,trunc(i/nim[1]))
       ii=i-ni*nim[1]+1
      x1=x[im[ii]]
      y1=y[im[ii]]
      z1=z[im[ii]]
     end 
     # ir=convert(Int64,trunc(nd[1]*rand()+1))  #...selects a random cell from the i_source selection done above
     #ijk = CartesianIndices(d)[i_source[ir]]
           ix=x1+rand()  #ijk[1]+rand()
           iy=y1+rand() #ijk[2]+rand()
           iz=z1+rand()  #]ijk[3]+rand()
           θ=pi*rand()
           φ=2*pi*rand()        
           vx=vc*sin(θ)*cos(φ)     #...random initial velocity vector (v=c)
           vy=vc*sin(θ)*sin(φ)
           vz=vc*cos(θ)
        p[1:3,ip].=[ix,iy,iz]
        p[4:6,ip].=[vx,vy,vz]
        p[7:9,ip].=[ix,iy,iz]

        if nE[1]==1 && E_initial[1]>0
        p[10,ip]=10^E_initial[1]  #...initial energy in eV 
        end 
        if nE[1]>1 && E_initial[1]>0
           ii=convert(Int64,round(nE[1]*rand()))
           if ii<1
           ii=1
           end 
           if ii>nE[1]
           ii=nE[1]
           end 
           p[10,ip]=10^(E_initial[ii])  #....this generates are random distribution of initial energies picked from the given energy bins.
           end 

           if nE[1]>1 && E_initial[1]<0
              Er=(E_initial[3])*rand()
              Ein_random=E_initial[2]+Er[1]
              p[10,ip]=10^(Ein_random)  #....this generates are random distribution of initial energies picked from the given energy bins.
        
        
        end 
      end 

   d=nothing
  return p
  end


 #...Boris pusher   
function integ_kdk(p1,dvx,pb,dt,qm,scale,ng,courant,iqm,energy,dEdt,zz)

    vp2=p1[4:6]
    modv=(p1[4]^2+p1[5]^2+p1[6]^2)  #...velocity module before kick (must be=c)
    vp2.+=dvx.*dt*0.5  #..kick   - velocity is updated for half a step 
   
   #...enforce energy conservation
   diss=(vp2[1]^2+vp2[2]^2+vp2[3]^2) 
   eps2=(modv-diss)/(modv) #..energy dissipation by the numerical scheme
   vp2.*=sqrt(modv/diss)   #..ensures the velocity module after kick is =c
 
   #...positions with enforced periodic BC (see below)
   p1[1]+=vp2[1]*dt/(scale)  #...drift
   p1[2]+=vp2[2]*dt/(scale)  #...drift
   p1[3]+=vp2[3]*dt/(scale)  #...drift

   #...positions without enforced periodic BC
   p1[7]+=vp2[1]*dt/(scale) #...drift
   p1[8]+=vp2[2]*dt/(scale) #...drift
   p1[9]+=vp2[3]*dt/(scale)  #...drift


      #...periodic boundary conditions
     if p1[1] < 1 
        p1[1]=ng+p1[1]-1
     end 
     if p1[2] < 1 
        p1[2]=ng+p1[2]-1
     end 
     if p1[3] < 1 
        p1[3]=ng+p1[3]-1
     end 
     if p1[1]>ng
     p1[1]=p1[1]-ng+1
    end  
    if p1[2]>ng
        p1[2]=p1[2]-ng+1
       end  
       if p1[3]>ng
        p1[3]=p1[3]-ng+1
       end  

       if isnan(p1[1])==1 || isnan(p1[2])==1 || isnan(p1[3])==1
        p1[1:10].=0.0
        end 

       i1=convert(Int64,trunc(p1[1]))
       i2=convert(Int64,trunc(p1[2]))
       i3=convert(Int64,trunc(p1[3]))

     dvx2=qm*cross_prod(vp2,pb)/vc
     if isnan(dvx2[1]) || isnan(dvx2[2])  || isnan(dvx2[3])   #...in case of some odd Nan numbers
       dvx2.=0.0 
     end   
     p1[4:6].=vp2+dvx2*dt*0.5   #...kick
   
   #...enforce energy conservation
   diss=(p1[4]^2+p1[5]^2+p1[6]^2)
   eps2=abs((modv-diss)/(modv)) #...energy dissipation by the numerical scheme
   p1[4:6].*=sqrt(modv/(diss))

    ebin=1
    neb=size(energy)

     #...the UHECR energy is redshifted to consider interaction with (1+z) more energetic photons, since we use z=0 EBL and CMB. 
     #...however if cosmo=0, we automatically get z=0 here 
     if p1[10]>maximum(energy)     #...in this case we don't have appropriate losses to evolve the particles 
      println("maximum energy exceeded")
      error()
      end 
     
   @inbounds for ie in eachindex(energy) #....locates the energy bin of the particle at this moment
   if p1[10]*(1+zz)>energy[ie] && p1[10]*(1+zz)<energy[ie+1]
   ebin=ie
   break 
   end 
end 

#....energy losses are computed at the end of the step, and energy is reduced accordingly 
#....loss terms include adiabatic expansion dEdt[ebin,1] and photon-pair and positron production with CMB+EBL dEdt[ebin,2], with z=0 data. 
#...for redshift dependences, see e.g. https://arxiv.org/pdf/hep-ph/0204357 eq.3-5 
  if p1[10]>0
dtloss=p1[10]/dEdt[ebin]
   ntime=convert(Int64,trunc(dt/dtloss))
   #....losses are computed in a single step, or in a loop, if dE/dt is large 
   Ez=sqrt(cOmegaM*(1+zz)^3+(1-cOmegaM))  #..adiabatic expansion term 
   if ntime<=1
      p1[10]=p1[10]*(1-((1+zz)^2*dEdt[ebin,2]+Ez*dEdt[ebin,1])*dt)  #...effect of energy loss from a particle with energy energy[ide]. Adiabatic losses are z-independent, while the losses with photon bacgrounds scale with (1+z)^2
   end 

   if ntime>=2
   @inbounds   for tt in 1:ntime
      p1[10]=p1[10]*(1-((1+zz)^2*dEdt[ebin,2]+(1+zz)^3*dEdt[ebin,1])*dtloss) 
   end 
end 
end 

    return p1[1:10]
      end
   
function move_CR(p::Array{Float64,2},t::Int64,skip_path::Int64,scale::Float64,courant::Float64,dt::Float64,dx::Float64,i1::Int64,i2::Int64,j1::Int64,j2::Int64,l1::Int64,l2::Int64,path::Array{Float64,3},cd::Float64,cv::Float64,cb::Float64,energy::Array{Float64},dEdt::Array{Float64},Z::Int64,zed::Float64,ngen::Int64,inj::Int64,bx::Array{Float64},by::Array{Float64},bz::Array{Float64})


  #....we apply or not cosmological correction factors related to z 
     if cosmo==0   #...the redshift is fixed to 0 everywhere 
     zed=0.0
     end 
    if cosmo==1    #...the comoving scale is made proper at the given z  
     scale=scale/(1+zed)
    end 

    #...selection of atomic mass number based on the nuclear charge number
    #...this may be set in a more general way for all particles at once, but it does not cost time and it allows us to assign different particles a different Z if we want 
    if Z==1  #proton 
    A=1
    end     
    if Z==2  #He
    A=4
    end    
    if Z==7  #Ni    
    A=14
    end    
    if Z==26  #Fe
    A=56
    end     
    np1=size(p)
 
    ng=i2-i1+1

   npt=convert(Int64,trunc(np1[2]/ngen))
   nev=convert(Int64,npt*inj)
   @inbounds for i in 1:nev   #...only evolve the UHECRs injected so far - the other remain idle 
     zz=zed  #...we assume the redshift as constant in this interval 

    i1=convert(Int64,trunc(p[1,i]))
    i2=convert(Int64,trunc(p[2,i]))
    i3=convert(Int64,trunc(p[3,i]))

    #...periodic boundary conditions are enforced here 
     if i1 < 1 
    i1=ng+i1
     end 
     if i2 < 1 
     i2=ng+i2
     end 
    if i3 < 1 
    i3=ng+i3
     end 
     if i1>ng
     i1=i1-ng
     end  
     if i2>ng
     i2=i2-ng
     end  
     if i3>ng
     i3=i3-ng
     end  
     vp=p[4:6,i]

       pb=[bx[i1,i2,i3],by[i1,i2,i3],bz[i1,i2,i3]].*(1+zed)^2   #.....magnetic field in the particle reference frame 
  
       γ=p[10,i]*evtoerg/(A*prest)     #....Lorentz factor of particles.

      qm=Z*qe/(A*mp*γ)  #...mass, charge and gamma of the particle to be set here 
      iqm=1/(qm)
      dvx=qm*cross_prod(vp,pb)#/vc   #....acceleration from Lorentz force

      rl=iqm/(sqrt(pb[1]^2+pb[2]^2+pb[3]^2))
      p1=p[1:10,i]

       pnew=integ_kdk(p1,dvx,pb,dt,qm,scale,ng,courant,iqm,energy,dEdt,zz)  #...integrator of particle motion with kick-drift-kick method in the Borish pusher
       p[1:10,i].=pnew
 
    #....we write in the path[] file (to be written on disk) only one step every skip_path, to save memory 
    it_path=convert(Int64,trunc(t/skip_path))
       if t/skip_path==it_path 
        path[i,1,it_path]=p[1,i]
         path[i,2,it_path]=p[2,i]
         path[i,3,it_path]=p[3,i]
         path[i,4,it_path]=p[10,i]
         path[i,5,it_path]=p[7,i]
         path[i,6,it_path]=p[8,i]
         path[i,7,it_path]=p[9,i]
         path[i,8,it_path]=sqrt(pb[1]^2+pb[2]^2+pb[3]^2)  #...physical magnetic field amplitude
         path[i,9,it_path]=zed  #...redshift
         path[i,10,it_path]=t*dt    #...time since the start of the simulation, in seconds
         
   end 


 end
 
return path
end


function losses(Z,main,z)   #....old function using tabulated data

    #....tabulated loss function 
    if Z==1 
        file_losses=string(main,"/proton_losses_all.dat")
       end 
       if Z==2
         file_losses=string(main,"/helium_losses_all.dat")
        end 
        if Z==7 
         file_losses=string(main,"/nitrogen_losses_all.dat")  #...file corrotto, devo sistemare 
        end  
        if Z==26 
         file_losses=string(main,"/iron_losses_all.dat")
        end 
          a=readdlm(file_losses)
          energy=a[:,1]
          loss_scale=a[:,2]
          dEdt=similar(energy)
          timesE=similar(energy)
          for e in eachindex(energy)
          dEdt[e]=vc/(loss_scale[e]*Mpctocm)     
          timesE[e]=1/dEdt[e]/3e16
          end 
        #..............................
        #...plot just to check the used loss function, can be omitted 
       lfs=14
       xfs=12
       xtfs=12
       plot(energy,timesE,color="blue",label=label=string("Z=",Z),line=:solid,dpi=500,linewidth=2.5,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
       title!(string("loss timescale, z=0"),fonts=20)
       yaxis!(L"[Gyr]",:log10,(0.001,1e2),fonts=20)
       xaxis!(L"E[eV]",:log10,(1e18,1e23),fonts=20)
       model=""
    
       filep1=string(main,"losses_timescale.png")
       savefig(filep1)
   
       plot(energy,dEdt,color="blue",label=string("Z=",Z),line=:solid,dpi=500,linewidth=2.5,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
       title!(string("loss scale, z=0"),fonts=20)
       yaxis!(L"dE/dt",:log10,fonts=20)
       xaxis!(L"E[eV]",:log10,(1e18,1e23),fonts=20)
       model=""
    
       filep1=string(main,"dEdt.png")
       savefig(filep1)
   
 return energy,dEdt
        end 



function lossesC(Z,main)  #...new function by C. Evoli, only for protons

   #....tabulated loss function 
   if Z==1 
       file_losses=string(main,"/SimProp_losses_proton.txt") 
      end 
     
         a=readdlm(file_losses)
         energy=a[:,1]
         nE=size(energy)
         loss_adiabatic=a[:,2]  #adiabatic
         loss_pair_CMB=a[:,3]  #pair production  CMB
         loss_pair_EBL=a[:,4]  #pair production - EBL
         loss_pp_CMB=a[:,5]  # photo pion production - CMB 
         loss_pp_EBL=a[:,6]  #photo pion production - EBL 
         dEdt=Array{Float64}(undef,nE[1],2)
         timesE=similar(energy)
         for e in eachindex(energy)
         dEdt[e,1]=(loss_adiabatic[e])/3e16
         dEdt[e,2]=(loss_pp_CMB[e]+loss_pp_EBL[e]+loss_pair_CMB[e]+loss_pair_EBL[e])/3e16
         timesE[e]=1/(dEdt[e,1]+dEdt[e,2])/3e16
         end 
       #..............................
       #...plot just to check the used loss function, can be omitted 
      lfs=14
      xfs=12
      xtfs=12
      plot(energy,timesE,color="blue",label=label=string("Z=",Z),line=:solid,dpi=500,linewidth=2.5,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
      title!(string("loss timescale, z=0"),fonts=20)
      yaxis!(L"Gyr",:log10,(1e-6,1e2),fonts=20)
      xaxis!(L"E[eV]",:log10,(1e14,1e21),fonts=20)
      model=""
   
      filep1=string(main,"losses_timescale_new.png")
      savefig(filep1)
  
      plot(energy,dEdt[:,1].+dEdt[:,2],color="blue",label=string("Z=",Z),line=:solid,dpi=500,linewidth=2.5,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
      title!(string("loss scale, z=0"),fonts=20)
      yaxis!(L"(dE/dt)",:log10,fonts=20)
      xaxis!(L"E[eV]",:log10,(1e14,1e21),fonts=20)
      model=""
   
      filep1=string(main,"dEdt_new.png")
      savefig(filep1)

return energy,dEdt
       end 

       function define_times(courant,scale,zfin,zin,ch,cOmegaM)

    #....DEFINITION OF THE FINAL ARRAY WITH ALL PARTICLE INFORMATION
    dt=courant*scale/vc
    max_it=convert(Int64,trunc(time_tot/dt))      #...maximum number of iterations
    println("going to evolve UHECR for ",max_it," iterations, to cover ",time_tot/Gyr, "Gyr of evolution") 

    times=Array{Float64}(undef,max_it) 
    zed=Array{Float64}(undef,max_it)
  
    #...DEFINITION OF TIMESTEPS
    #...we use the approximated inversion between time and scale factor as in https://academic.oup.com/mnras/article/505/2/2764/6289943, eq.15-17 
    dz=(zin-zfin)/max_it
#    time_tilde=0.666/(ch*100*sqrt(1-cOmegaM))/Gyr
    time_tilde=(6.519/ch)*(1/sqrt(1-cOmegaM))
    @inbounds for i in 1:max_it #...we define an appropriate array of times and redshifts 
    times[i]=dt*i/Gyr
    a=(sqrt(cOmegaM/(1-cOmegaM))*sinh((Δt/Gyr+times[i])/time_tilde))^0.6666         
    zed[i]=1/a-1     
   end 

  return times,zed,dt,max_it 

end 


function inject_new_UHECR(root0,root_halos,snapiz,p,i1,i2,j1,j2,l1,l2,mass_source,E_initial,np,inj,ngen,n)
       
    file1=string(root0,"minus1.0_08RD00",snapiz,".cpu0000")   #...ENZO snapshots with 3D data
    filecat=string(root_halos,snapiz,"_minus1.0_08_halof_200_new.dat")   #...halo catalogs 

   # p=assign_CR_dens_z(p,i1,i2,j1,j2,l1,l2,file1,dsource,E_initial,np,inj,ngen)   #...assign UHECR based on density 
    p=assign_CR_halo_z(p,i1,i2,j1,j2,l1,l2,filecat,mass_source,E_initial,np,inj,ngen,n)   #...assign UHECR based on density  
    inj0=inj
    ntr=size(p)    #....number of UHECR - it depends on the dsource density choosen above 
    np=ntr[2]
    npt=convert(Int64,trunc(np/ngen))   #....number of new UHECR to be injected at this step 
    nj1=(inj-1)*npt+1
    nj2=nj1+npt-1

   #....we plot the initial location of UHECR
    plo=plot(p[1,nj1:nj2],p[2,nj1:nj2],seriestype=:scatter,ms=1,label="",grid=false,aspect_ratio=1.0)  #...plot initial location of CRs
    xlims!(i1, i2)
    ylims!(j1,j2)
    filep1=string(root_out,"_initial_map_newU",snapiz,"_cosmo",cosmo,".png")
        savefig(filep1)

    if use_syntheticB==1   #...using a digitally-generated 3D B-field 
      #.............this is a temporary hardcoding for some testing, but it is ugly and should be removed 
      fold_ic="/Users/francovazza/Library/CloudStorage/Dropbox/Julia_prog/B_init/150Mpc/"
      filex=string(fold_ic,"Bx_-1.0")
      filey=string(fold_ic,"By_-1.0")
      filez=string(fold_ic,"Bz_-1.0")
      bx=h5read(filex,"Bx",(i1:i2,j1:j2,l1:l2)) 
      by=h5read(filey,"By",(i1:i2,j1:j2,l1:l2)) 
      bz=h5read(filez,"Bz",(i1:i2,j1:j2,l1:l2)) 
      zin=30.0
     global  bx./=(1+zin)^2    #...comoving B-field in Gauss , zin depends on the specific external box we are reading here 
     global  by./=(1+zin)^2
     global  bz./=(1+zin)^2
    bx=convert(Array{Float64,3},bx.*0.0)
    by=convert(Array{Float64,3},by.*0.0)
    bz=convert(Array{Float64,3},bz.*0.0)
      end 

    if use_syntheticB==0    #...default use with B-fields evolved by ENZO 
      bx=h5read(file1,"Grid00000001/Bx",(i1:i2,j1:j2,l1:l2))
      by=h5read(file1,"Grid00000001/By",(i1:i2,j1:j2,l1:l2))
      bz=h5read(file1,"Grid00000001/Bz",(i1:i2,j1:j2,l1:l2))     
       global bx*=cb   #...this field is comoving. it is turned into physical, x(1+z)^2, inside the propagation routine
       global by*=cb
       global bz*=cb
       println(typeof(bx))
      end 


      return p,bx,by,bz
end 




function plot_map_WBC(dx,np,path,E_initial,Z,root_out,cosmo,tag)

          #1) static plot without boundary conditions
          xs=dx*1e-3  #...to convert positions from cells to Mpc
          jump=10
          #the following plot can use a lot of memory, do for i in 1:jump:np with jump = any integer number to plot only 1 every jump particles
          @inbounds for i in 1:jump:np-1   #....all particle trajectories are plotted without periodic BC 
          if i==1 
          plo=plot(path[i,5,:]*xs,path[i,6,:]*xs,label="",aspect_ratio = :equal ,dpi=1000,lw=0.5,grid=false)
          end 
          io=findall((path[i,5,:].>0) .| (path[i,5,:].<0))
#          io=findall(abs.(path[i,5,:]). > 0) 
          plot!(path[i,5,io]*xs,path[i,6,io]*xs,label="",lw=0.5,alpha=0.5)#,seriestype=:scatter,ms=0.4)
          end
          tit=string(" E0=",E_initial)
          if E_initial[1]==-0
            tit=L" N(E) \propto E^{-1}"
          end 
          title!(string("Z=",Z," ,",tit," eV"),fonts=20)
              yaxis!("[Mpc]",fonts=20)
              xaxis!("[Mpc]",fonts=20)
              filep1=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_noBC_cosmo",cosmo,"_",tag,".png")
              savefig(filep1)

        end 



    function make_gif(skip_path,npath,E_initial,path,nz,Z,root_out,n,cosmo,tag)
          #2) animated gif with periodic boundary conditions
          #...plotting evolution of UHECR in the periodic BC as an animated gif 
        it=0
        xs=dx*1e-3  #...to convert positions from cells to Mpc
   
          #...this plotting function in general produces way too many frames and it wastes some time
          #   nframe=npath  if one wants to produce the full movie with all steps in the path file
          #   nframe<npath  if one wants to limit the gif to a smaller number of frames, like npath=100 for the first 100 steps only  
            nframe=npath
            anim=@animate for t in 1:skip_path:nframe*skip_path   #....plotting positions 
               it+=1
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
          plot!(path[g1:g2,1,it]*xs,path[g1:g2,2,it]*xs,label="",line=0,lw=0.5,marker=:circle,ms=1.2,seriestype=:scatter,color=cg[gg],alpha=0.5)
        end 
        title!(string("Z=",Z," ,",tit," eV, time= ",trunc(times[it*skip_path],digits=4),"Gyr"),fonts=20)
        yaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        xaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
        end 
        file_gif=string(root_out,"_UHECR_path_color_map_",E_initial[1],"_Z_",Z,"_BC_redshift_newU_cosmo",cosmo,"_",tag,".gif")
        gif(anim,file_gif,fps=10)    
        
    end 


    function plot_spectrum(path,np,root_out,E_initial,Z,cosmo,tag)


         #....plotting the final spectrum 
         nspec=10
         mie=16.0   #...logE
         mae=21.0
         bie=(mae-mie)/nspec
         xe=Array{Float64}(undef,nspec)
         spec=Array{Float64}(undef,nspec,npath)
         spec.=0.0
         iobs=npath  #...snapshot where we produce the final spectrum 
          xe.=0.0
           @inbounds  for i in 1:nspec
            xe[i]=10.0^(mie+bie*(i-0.5))
            end

      @inbounds for j in 1:npath
      @inbounds @simd for i in 1:np 
      en=path[i,4,j]
      if en>10.0^mie && en<10.0^mae  
      ie=convert(Int64,trunc((log10(en)-mie)/bie))+1
     spec[ie,j]+=1.0
      end 
      end 
      end 


             e1=17.0
             e2=21.0
            plot(xe,spec[:,1],line=:solid,dpi=300,lw=1,grid=false,alpha=0.5,label=string("z=",trunc(path[1,9,1],digits=6)),color="black")
             plot!(xe,spec[:,iobs],line=:solid,lw=2,alpha=0.5,label=string("z=0"),color="red")
             yaxis!("N(E)",:log10,(1,np),fonts=20)
             xaxis!(L"log_{10}E[eV]",:log10,(10^e1,3*10^e2),fonts=20,xticks=[1e17,1e18,1e19,1e20,1e21])
         
            file=string(root_out,"_UHECR_spectum_",E_initial[1],"_Z_",Z,"_newU_cosmo",cosmo,"_",tag,".png")
            savefig(file)


    end 



    function plot_spectrum_suppression(path,np,root_out,E_initial,Z,cosmo,tag,root_halos)


        #....plotting the final spectrum 
        nspec=10
        mie=16.0   #...logE
        mae=21.0
        bie=(mae-mie)/nspec
        xe=Array{Float64}(undef,nspec)
        spec=Array{Float64}(undef,nspec,npath)
        spec.=0.0
        iobs=npath  #...snapshot where we produce the final spectrum 
         xe.=0.0
          @inbounds  for i in 1:nspec
           xe[i]=10.0^(mie+bie*(i-0.5))
           end

     @inbounds for j in 1:npath
     @inbounds @simd for i in 1:np 
     en=path[i,4,j]
     if en>10.0^mie && en<10.0^mae  
     ie=convert(Int64,trunc((log10(en)-mie)/bie))+1
    spec[ie,j]+=1.0
     end 
     end 
     end 

         spec_obs=Array{Float64}(undef,nspec,npath,2)
         spec_obs.=0.0
         mass_obs=7e11      
         mass_source=1e13   #...set this to what is used in the code  
         #....select observers at z=0
           filecat=string(root_halos,"21_minus1.0_08_halof_200_new.dat")   #...halo catalogs 
           a=readdlm(filecat)
           imass=findall((a[:,5].>=mass_obs) .& (a[:,5].<=mass_source))# masobs_source &&. a[:,5].>=1e12)
           println("number of observers=",size(imass))
           mass=a[imass,5]
           radius=a[imass,4]
           x=a[imass,1]
           y=a[imass,2]
           z=a[imass,3]
           x.*=n
           y.*=n
           z.*=n
           x.+=(-i1+1)
           y.+=(-j1+1)
           z.+=(-l1+1)      
           im=sortperm(mass,rev=true)    
           nim=size(im)
           radius=3  #...observer radius in unit CELLS (not Mpc)
           iobs=npath-1 #...snapshot where the observer is placed. better not the very last snapshot to avoid UHECRs injected in the last timestep 
           @inbounds for ii in eachindex(im)#...loop over observers 
             ox1=x[im[ii]]
             oy1=y[im[ii]]
             oz1=z[im[ii]]

             @inbounds for j in 1:np    #...loop over particles s
              dist=sqrt( (x1[j,iobs]-ox1)^2+(y1[j,iobs]-oy1)^2+(z1[j,iobs]-oz1)^2)
               if dist < radius 
                 en=E1[j,iobs]
                 if en>10.0^mie && en<10.0^mae  
                ie=convert(Int64,trunc((log10(en)-mie)/bie))+1
                  spec_obs[ie,1,1]+=1.0
                 end 
              end 

              dist=sqrt( (x2[j,iobs]-ox1)^2+(y2[j,iobs]-oy1)^2+(z2[j,iobs]-oz1)^2)
              if dist < radius #&& dist>1
                en=E2[j,iobs]
                if en>10.0^mie && en<10.0^mae  
               ie=convert(Int64,trunc((log10(en)-mie)/bie))+1
                 spec_obs[ie,1,2]+=1.0
                end 
             end 
             end 


           end 

      e1=17.0
      e2=21.0
      E_ref=19.0    #.....reference log(E) at which we normalise the total and the observed spectrum to have the same value

      Eref=convert(Int64,trunc((E_ref-mie)/bie))
      norm0=sum(spec[Eref,iobs,1])/sum(spec[Eref,iobs,2])
      norm1=sum(spec[Eref,iobs,1])/sum(spec_obs[Eref,1,1])
      norm2=sum(spec[Eref,iobs,1])/sum(spec_obs[Eref,1,2])
      
plot(xe,spec[:,iobs,1],line=:solid,dpi=300,lw=1,grid=false,alpha=0.5,label=string("B"),color="black")
#plot!(xe,norm0*spec[:,iobs,2],line=:dash,lw=2,alpha=0.5,label=string("B=0"),color="black")
plot!(xe,norm1*spec_obs[:,1,1],line=:solid,lw=1,alpha=0.5,label=string("obs, B"),color="red")
#plot!(xe,norm2*spec_obs[:,1,2],line=:dash,lw=1,alpha=0.5,label=string("obs, B=0"),color="red")
  

                 yaxis!("N(E)",:log10,(1e2,np),yticks=[1,10,1e2,1e3,1e4,1e5],fonts=20)
                 xaxis!(L"log_{10}E[eV]",:log10,(3*10^e1,3*10^e2),fonts=20,xticks=[1e17,1e18,1e19,1e20,1e21])
          
            file=string(root_out,"_UHECR_spectum_cfr_",tag,".png") 
            savefig(file)

          sup=similar(spec_obs)
          sup.=0.0

          for i in 1:nspec 
           sup[i,1,1]=(norm1*spec_obs[i,1,1])/spec[i,iobs,1]
          end 

                     plot(xe,sup[:,1,1],line=:solid,dpi=300,lw=1,grid=false,alpha=0.5,color="black",label="z=0")
                     yaxis!(L"G(E)=N(E)_{obs}/N(E)",fonts=20,(0,1.2))
                     xaxis!(L"log_{10}E[eV]",:log10,(3*10^e1,3*10^e2),fonts=20,xticks=[1e17,1e18,1e19,1e20,1e21])
                   
                      file=string(root_out,"_UHECR_spectum_suppressions_",tag,".png")
                      savefig(file)
            
   end 






     function write_path(root_out,E_initial,Z,cosmo,tag,path)

      #...for safety, HDF5 files must be manually deleted before overwriting - otherwise the code stops here
      filep1=string(root_out,"path_",E_initial[1],"Z",Z,"_spec_redshift_cosmo",cosmo,"_",tag,".hdf5")
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


     end 