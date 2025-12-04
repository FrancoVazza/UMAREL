  
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
    d=  h5read(file1,"Density",(i1:i2,j1:j2,l1:l2)) 
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
               Ein_random=17+Er[1]
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
              Er=E_initial[3]*rand()
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
   @inbounds for ie in eachindex(energy) #....locates the energy bin of the particle at this moment
   if p1[10]>energy[ie] && p1[10]<energy[ie+1]
   ebin=ie
  
   break 
end 
end 


#....energy losses are computed at the end of the step, and energy is reduced accordingly 
   dtloss=p1[10]/dEdt[ebin]
   ntime=convert(Int64,trunc(dt/dtloss))
   #....losses are computed in a single step, or in a loop, if dE/dt is large 
   if ntime<=1
      p1[10]=p1[10]*(1-((1+zz)^2*dEdt[ebin,2]+dEdt[ebin,1])*dt)  #...effect of energy loss from a particle with energy energy[ide]. Adiabatic losses are z-independent, while the losses with photon bacgrounds scale with (1+z)^2
   end 

   if ntime>=2
   @inbounds   for tt in 1:ntime
      p1[10]=p1[10]*(1-((1+zz)^2*dEdt[ebin,2]+dEdt[ebin,1])*dtloss) 
#      p1[10]=p1[10]*(1-dEdt[ebin]*dtloss)
   end 
end 

    return p1[1:10]
      end
   
function move_CR(p::Array{Float64,2},t::Int64,skip_path::Int64,scale::Float64,courant::Float64,dt::Float64,dx::Float64,i1::Int64,i2::Int64,j1::Int64,j2::Int64,l1::Int64,l2::Int64,path::Array{Float64,3},cd::Float64,cv::Float64,cb::Float64,energy::Array{Float64},dEdt::Array{Float64},Z::Int64,zed::Float64,ngen::Int64,inj::Int64,file1,bx::Array{Float64},by::Array{Float64},bz::Array{Float64})
#...selection of atomic mass number based on the nuclear charge number
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

       pb=[bx[i1,i2,i3]*(1+zz)^2,by[i1,i2,i3]*(1+zz)^2,bz[i1,i2,i3]*(1+zz)^2]   #.....physical B-field (1+z)^2
  
       γ=p[10,i]*evtoerg/(A*prest)     #....Lorentz factor of particls. Notice it should be changed for UHECR with a higher composition!!!!!! 

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


function losses(Z,main,z)

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



function lossesC(Z,main)

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