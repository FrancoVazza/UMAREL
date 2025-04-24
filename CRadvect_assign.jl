  
function cross_prod(a,b)
   c=similar(a)
   c[1] = a[2]*b[3]-a[3]*b[2]
   c[2] = a[3]*b[1]-a[1]*b[3]
   c[3] = a[1]*b[2]-a[2]*b[1]
   return c 
end

  
  function assign_CR_dens(i1,i2,j1,j2,l1,l2,file1,dsource,E_initial)
    #...this function assigns the initial position of UHECRs within cells denser than something ("halos")
    #...depending on the dsource density threshold, more or less UHECRs are injected 
   d=h5read(file1,"Density",(i1:i2,j1:j2,l1:l2))   
   dCR=dsource*mean(d)  #...to be changed with absolute density threshold? 
   d0=similar(d)
   println("density threshold=", dCR)
@inbounds   for i in eachindex(d)
   d0[i]=trunc(d[i]/dCR)
   end  
   id=findall((d0.>= 1))
   totd=sum(d0[id])
    
   np=convert(Int64,trunc(totd))
   println(np)

   p=Array{Float64}(undef,10,np)
   p.=0.0

   npart=0
   Random.seed!(123)
     @inbounds for l in 1:n
      @inbounds for j in 1:n
         @inbounds for i in 1:n
            nt=convert(Int64,trunc(d[i,j,l]/dCR))#   number of UHECR to be generated in the cell
            if nt >=1
               npart+=nt
            @inbounds @simd for c in 1:nt #....loops on sorted cels 
            ix=i+rand()
            iy=j+rand()
            iz=l+rand()
            θ=pi*rand()
            φ=2*pi*rand()        
            vx=vc*sin(θ)*cos(φ)     #...random initial velocity vector (v=c)
            vy=vc*sin(θ)*sin(φ)
            vz=vc*cos(θ)
         p[1:3,npart-c+1].=[ix,iy,iz]
         p[4:6,npart-c+1].=[vx,vy,vz]
         p[7:9,npart-c+1].=[ix,iy,iz]

         if nE[1]==1
         p[10,npart-c+1]=10^E_initial[1]  #...initial energy in eV 
         end 
         if nE[1]>1
            ii=convert(Int64,trunc(nE[1]*rand()+1))
         p[10,npart-c+1]=10^(E_initial[ii])  #....this generates are random distribution of initial energies from 10^18 to 10^(18+E_initial)
         end 
         end 
         end 
      end
   end
end
    
    d=nothing
   return p
   end
 

 #...Boris pusher   
function integ_kdk(p1,dvx,pb,dt,qm,scale,ng,courant,iqm,energy,dEdt)
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
   p1[7]+=vp2[1]*dt/(scale)  #...drift
   p1[8]+=vp2[2]*dt/(scale)  #...drift
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
      p1[10]=p1[10]*(1-dEdt[ebin]*dt)  #...effect of energy loss from a particle with energy energy[ide]
   end 

   if ntime>=2
   @inbounds   for tt in 1:ntime
      p1[10]=p1[10]*(1-dEdt[ebin]*dtloss)
   end 
end 

    return p1[1:10]
      end
   
function move_CR(p::Array{Float64,2},max_it::Int64,skip_path::Int64,scale::Float64,courant::Float64,dt::Float64,dx::Float64,i1::Int64,i2::Int64,j1::Int64,j2::Int64,l1::Int64,l2::Int64,path::Array{Float64,3},cd::Float64,cv::Float64,cb::Float64,energy::Array{Float64},dEdt::Array{Float64},Z::Int64)
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
    np=size(p)
    bx=h5read(file1,"Bx",(i1:i2,j1:j2,l1:l2))
    by=h5read(file1,"By",(i1:i2,j1:j2,l1:l2))
    bz=h5read(file1,"Bz",(i1:i2,j1:j2,l1:l2))
    
    ng=i2-i1+1

    bx*=cb
    by*=cb
    bz*=cb

#....fixed magnetic field for testing. For a negligible B-field, set bx.=1e-20 etc,, for bx.=1e-9 for an important one 
#    bx.=2e-9
#    by.=2e-9
#    bz.=2e-9
    

 println("total time of propagation=",time_tot)

@inbounds for i in 1:np[2]   #...loop over all particles 
@inbounds for it in 1:max_it    #....time integration for each particle 

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

    pb=[bx[i1,i2,i3],by[i1,i2,i3],bz[i1,i2,i3]]
    γ=p[10,i]*evtoerg/(A*prest)     #....Lorentz factor of particls. Notice it should be changed for UHECR with a higher composition!!!!!! 

      qm=Z*qe/(A*mp*γ)  #...mass, charge and gamma of the particle to be set here 
    iqm=1/(qm)
    dvx=qm*cross_prod(vp,pb)#/vc   #....acceleration from Lorentz force
  
   rl=iqm/(sqrt(pb[1]^2+pb[2]^2+pb[3]^2))
   p1=p[1:10,i]

    pnew=integ_kdk(p1,dvx,pb,dt,qm,scale,ng,courant,iqm,energy,dEdt)  #...integrator of particle motion with kick-drift-kick method in the Borish pusher
    p[1:10,i].=pnew
 
    #....we write in the path[] file (to be written on disk) only one step every skip_path, to save memory 
    it_path=convert(Int64,trunc(it/skip_path))
    if it/skip_path==it_path 
    path[i,1,it_path]=p[1,i]
    path[i,2,it_path]=p[2,i]
    path[i,3,it_path]=p[3,i]
    path[i,4,it_path]=p[10,i]
    path[i,5,it_path]=p[7,i]
    path[i,6,it_path]=p[8,i]
    path[i,7,it_path]=p[9,i]
    end 


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
          times=similar(energy)
          for e in eachindex(energy)
          dEdt[e]=vc/(loss_scale[e]*Mpctocm)     
          times[e]=1/dEdt[e]/3e16
          end 
        #..............................
        #...plot just to check the used loss function, can be omitted 
       lfs=14
       xfs=12
       xtfs=12
       plot(energy,times,color="blue",label=label=string("Z=",Z),line=:solid,dpi=500,linewidth=2.5,alpha=0.7,grid=false,legendfontsize=lfs,yguidefontsize=xfs,xguidefontsize=xfs,xtickfontsize=xtfs,ytickfontsize=xtfs)
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