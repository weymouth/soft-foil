using WaterLily,Plots,Roots,StaticArrays

function plot_vorticity(sim)
    @inside sim.flow.σ[I] = sum(WaterLily.ϕ(i,CartesianIndex(I,i),sim.flow.μ₀) for i ∈ 1:2)/2
    @inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u)*sim.flow.σ[I]*sim.L / sim.U
    contourf(sim.flow.σ',dpi=300,
			 color=palette(:BuGn), clims=(-10, 10), linewidth=0,
			 aspect_ratio=:equal, legend=false, border=:none)
end

function foil(L;d=0,Re=38e3,U=1)
	# Define symmetric foil shape with cubic Bezier
    P0,P1,P2,P3 = SA[0,0],SA[0,0.1],SA[.5,0.12],SA[1.,0.]
    A,B,C,D = P0,-3P0+3P1,3P0-6P1+3P2,-P0+3P1-3P2+P3
    curve(t) = A+B*t+C*t^2+D*t^3
    dcurve(t) = B+2C*t+3D*t^2

	# Find distance to the symmetric foil
    candidates(X) = union(0,1,find_zeros(t -> (X-curve(t))'*dcurve(t),0,1,naive=true,no_pts=3,xatol=0.01))
    function distance(X,t)
        V = X-curve(t)
        copysign(√sum(abs2,V),[V[2],-V[1]]'*dcurve(t))
    end
    distance(X) = L*argmin(abs, distance(X,t) for t in candidates(X))
    
    # Define pitching motion mapping
    ω = 2π*U/L # reduced frequency k=π
    center = SA[1,1] # foil placement in domain
    pivot = SA[.1,0] # pitch location from center
    function map(x,t)
        α = 6π/180*cos(ω*t)
        SA[cos(α) sin(α); -sin(α) cos(α)] * (x/L.-center-pivot) .+ pivot
    end

    # Define sdf to deflected foil and Simulation
    deflect(x) = max(0,x-0.3)^2/0.7^2
    sdf(x,t) = distance(SA[x[1],abs(x[2]+d*deflect(x[1]))])
    Simulation((4L+2,2L+2),[U,0],L;ν=U*L/Re,body=AutoBody(sdf,map))
end

sim = foil(256,d=0.12);
sim_step!(sim,4,remeasure=true);
plot_vorticity(sim)
savefig("flow_d_12.png")

@gif for t ∈ 4:1/24:(6-1/24)
    sim_step!(sim,t,remeasure=true)
    plot_vorticity(sim)
end