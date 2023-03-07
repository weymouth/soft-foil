using ProperOrthogonalDecomposition,HDF5,Printf
include("simulation.jl")
export foil,get_omega!,plot_vorticity

filename(L,d,tail) = @sprintf ".\\flow\\L%i_d%02i_%s" L 100d tail

function get_data(L,d;warmup=6::Int,collection=2::Int,dt=1/24,m=3)
    # Define the simulation
    println("Starting simulation:",filename(L,d,"dat"))
    sim = foil(L;d,m)

    # Run the warmup without saving
    for t in 1:warmup
        sim_step!(sim,t,remeasure=true);
        println("warmup time:",t," out of ",warmup)  
    end

    # Save a plot 
    get_omega!(sim)
    plot_vorticity(sim.flow.σ,limit=10)
    savefig(filename(L,d,"dat.png"))

    # Run data collection phase
    data = mapreduce(hcat,(warmup+dt):dt:(warmup+collection)) do t
        sim_step!(sim,t,remeasure=true)
        println("collection time:",round(t,digits=2)," out of ",warmup+collection)
        get_omega!(sim)
        reshape(sim.flow.σ,(:))
    end |> (x -> reshape(x,(:,m*L+2,Int(collection÷dt))))
    h5write(filename(L,d,"dat.h5"), "group\\data", data)
end

function get_POD(L,d;window=2L:5L,fname=filename(L,d,"pod.h5"))
    # Get data from file
    data = h5read(filename(L,d,"dat.h5"), "group\\data")
    n = size(data)

    # Take POD of windowed data
    data = reshape(data[window,:,:],(:,n[3]))
    res, singularvals = PODsvd!(data)
    data = reshape(res.modes,(:,n[2],n[3]))

    # (Over)Write to file
    h5open(fname, "w") do fid
        fid["group\\singularvals"] = singularvals
        fid["group\\coefficients"] = res.coefficients
        fid["group\\modes"] = data
    end
end

get_data(32,0,warmup=15,collection=10)