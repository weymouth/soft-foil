using WaterLily,Plots,Roots,StaticArrays,ProperOrthogonalDecomposition,HDF5,Printf

function get_omega!(sim)
    body(I) = sum(WaterLily.ϕ(i,CartesianIndex(I,i),sim.flow.μ₀) for i ∈ 1:2)/2
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u) * body(I) * sim.L / sim.U
end

plot_vorticity(ω; level=maximum(abs,ω)) =contourf(ω',dpi=300,
    color=palette(:BuGn), clims=(-level, level), linewidth=0,
    aspect_ratio=:equal, legend=false, border=:none)

function foil(L;d=0,Re=38e3,U=1,n=6,m=3)
	# Define symmetric foil shape with cubic Bezier
    P0,P1,P2,P3 = SA[0,0],SA[0,0.1],SA[.5,0.12],SA[1.,0.]
    A,B,C,D = P0,-3P0+3P1,3P0-6P1+3P2,-P0+3P1-3P2+P3
    curve(t) = A+B*t+C*t^2+D*t^3
    dcurve(t) = B+2C*t+3D*t^2

	# Signed distance to the symmetric foil
    candidates(X) = union(0,1,find_zeros(t -> (X-curve(t))'*dcurve(t),0,1,naive=true,no_pts=3,xatol=0.01))
    function distance(X,t)
        V = X-curve(t)
        copysign(√(V'*V),[V[2],-V[1]]'*dcurve(t))
    end
    distance(X) = L*argmin(abs, distance(X,t) for t in candidates(X))
    
    # Pitching motion around the pivot
    ω = 2π*U/L # reduced frequency k=π
    center = SA[1,m/2] # foil placement in domain
    pivot = SA[.1,0] # pitch location from center
    function map(x,t)
        α = 6π/180*cos(ω*t)
        SA[cos(α) sin(α); -sin(α) cos(α)]*(x/L-center-pivot) + pivot
    end

    # Define sdf to deflected foil and Simulation
    deflect(x) = max(0,x-0.3)^2/0.7^2
    sdf(x,t) = distance(SA[x[1],abs(x[2]+d*deflect(x[1]))])
    Simulation((n*L+2,m*L+2),[U,0],L;ν=U*L/Re,body=AutoBody(sdf,map))
end

filename(L,d,tail) = @sprintf "flow\\L%i_d%02i_%s" L 100d tail

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
    plot_vorticity(sim.flow.σ,level=10)
    savefig(filename(L,d,"dat.png"))

    # Run data collection phase
    data = mapreduce(hcat,(warmup+dt):dt:(warmup+collection)) do t
        sim_step!(sim,t,remeasure=true)
        println("collection time:",round(t,digits=2)," out of ",warmup+collection)
        get_omega!(sim)
        reshape(sim.flow.σ,(:))
    end |> (x -> reshape(x,(:,m*L+2,Int(collection÷dt))))
    h5write(filename(L,d,"dat.h5"), "group\\data", data)
    plot_vorticity(data[:,:,end],level=10)
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

    # Return a plot of the first mode as a sanity check
    plot_vorticity(data[:,:,1])
end