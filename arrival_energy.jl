using HDF5
using Plots
using Statistics
using LaTeXStrings
using Distributions
using StatsBase


#PESCO DALLA SEQUENZA DEGLI SNAPSHOP GLI ISTANTI DI TEMPO A CUI TUTTE LE PARTICELLE RAGGIUNGONO LA STESSA ENERGIA
#se io sono l'osservatore e mi vedo arrivare una particella con energia 1e20, voglio trovare quanta distanza ha percorso. (?) 
#1) Quindi devo selezionare dal path.hdf5 le celle in cui trovo l'energia 1e20, così ho il mio final timestep
#2) uso questo time step per calcolare la distanza percorsa (decine di Mpc)
      root0=string("/Users/francovazza/Desktop/data/DATA/RADGAL/LEONARDO/42.5Mpc/snap/")
       root_out=string(root0,"/out/")

root_out=string(root0,"/out/")
plot_dir = root0#joinpath(root0,"newplots")
#CREO FOLDER PER I PLOT 
if !isdir("newplots")
    mkdir("newplots")
end


main="/Users/francovazza/Dropbox/Julia_prog/UHECRpropa/MAGHOR/MAGHOR/" #..main folder containing ROGER functions
#...modules 
include(string(main,"/constants.jl"))
include(string(main,"CRadvect_assign.jl"))   #...external module with all relevant functions used for the transport of CRs
 
################################################################################################################

filep1 = string(root_out,"path_-1Z1_spec.hdf5")
  
px = h5read(filep1,"px")  #...X position without periodic BC 
py = h5read(filep1,"py")  #...Y position without periodic BC 
pz = h5read(filep1,"pz")  #...Z position without periodic BC 
energy = h5read(filep1,"E[eV]")  #...energy 
t = h5read(filep1,"time")         #...list of simulated timesteps
println("HO LETTO")
println(maximum(t))

nE=size(energy)
println(nE)
plot(t[:],energy[1,:],dpi=500,label="",title=string(nE[2]," particles"))
for i in 1:nE[2]
    plot!(t,energy[i,:],label="")
end 
    #yaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
 #       xaxis!("[Mpc]",(0,dx*n*1e-3),fonts=20)
yaxis!("[Mpc]",:log10,(5e17,1e22))
xaxis!("[Gyr]")

savefig(joinpath(plot_dir, "evolve_energy.png"))

t0=time()
p = plot()
for n in [18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5]

    println(".....................reaching",n," eV...........................")
    initial_cr, time_step = size(energy)

    distances = Array{Float64}(undef,initial_cr, time_step)
    distances.=0.0
 
    for i in 1: initial_cr

      #  println("------------------------------------------")
       # println("particle number ",i," has initial energy ",  energy[i, 1])
    
    
 #       id_final = findfirst(x -> isapprox(x, 10.0^ n, rtol=0.5), energy[i, :]) #la particella i raggiunge 1e20 al timestep id_final
        id_final=(findall(x->isapprox(x,10.0^n,rtol=0.25),energy[i,:]))
        
        ##################################################################
        if isnothing(id_final)
            println("particle",i,"does not reach",10^n)
            continue
        end


       # println("particle" ,i,"reaches", 10^n,"eV at ",id_final)

        ##################################################################

        #calcolo la distanza percorsa dal punto di iniezione iniziale al punto che corrisponde al timestep a cui raggiunge la data energia
        for j in eachindex(id_final)
            dx = px[i, j] - px[i, 1]
            dy = py[i, j] - py[i, 1]
            dz = pz[i, j] - pz[i, 1]

            distance =  sqrt(dx^2 + dy^2 + dz^2) * 0.0415 #in Mpc
          # println(distance)
            if distance>1e5 || isnan(distance)==1
                distance=0
                end  
            distances[i, j] = distance
        
         end
    
    
    #println("Particle $i travelled $distance Mpc at time step ", id_final)
    
    end
     
    println(maximum(distances))
   

    abs_max_d = convert(Int64, trunc(maximum(distances)))    

    println("the absolute maximum distances is", abs_max_d," Mpc")
    
    
    dist = Float64[]
    frac_cr = Float64[]

    @inbounds for d in 1:abs_max_d

        num_cr= count(i -> any(distances[i, :] .> d), 1 : initial_cr)


        append!(dist, d)
        append!(frac_cr, num_cr/initial_cr)


    end


    println("plotting...")

    plot!((dist), frac_cr, label = string(10^n," eV"), title = "Arrival energy", xlim =(0, 4), ylim = (0,1), xlabel = "log(distance [Mpc])", ylabel = "fraction of CR", theme =:rose_pine_dawn, legend =:outerright ) 
    xaxis!(:log10,(0.1,100))
end

println("ELAPSED TIME=",time()-t0)
savefig(joinpath(plot_dir, "arrival_energy.png"))

