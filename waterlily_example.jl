using WaterLily,Plots,Roots,StaticArrays

function foil(L;Re=38e3)
	# Define Bezier foil shape
    P0,P1,P2,P3 = SVector(0,0),SVector(0,0.1),SVector(0.5,0.12),SVector(1.,0.)
    A,B,C,D = P0,-3P0+3P1,3P0-6P1+3P2,-P0+3P1-3P2+P3
    curve(t) = A+B*t+C*t^2+D*t^3
    dcurve(t) = B+2C*t+3D*t^2

	# Define distance to foil
    candidates(X) = union(0,1,find_zeros(t -> (X-curve(t))'*dcurve(t),0,1,naive=true,no_pts=3,xatol=0.01))
    function distance(X,t)
        V = X-curve(t)
        return copysign(√sum(abs2,V),[V[2],-V[1]]'*dcurve(t))
    end
    distance(X) = argmin(abs,distance(X,t) for t in candidates(X))

	# Define body and simulation
	body = AutoBody((x,t) -> distance(SVector(x[1]/L-1,abs(x[2]/L-0.5))))
    return Simulation((4L+2,L+2),[1,0],L;ν=L/Re,body)
end

function plot_vorticity(sim)
	@inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
	contourf(sim.flow.σ',
			 color=palette(:BuGn), clims=(-5, 5), linewidth=0,
			 aspect_ratio=:equal, legend=false, border=:none)
end

sim = foil(64);
sim_step!(sim,10);
plot_vorticity(sim)
savefig("flow.png")