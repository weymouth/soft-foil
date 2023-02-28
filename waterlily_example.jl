using WaterLily,Plots
using LinearAlgebra: norm2

function circle(radius=8;Re=250,n=10,m=6)
    center, ν = radius*m/2, radius/Re
    body = AutoBody((x,t)->norm2(x .- center) - radius)
    Simulation((n*radius+2,m*radius+2), [1.,0.], radius; ν, body)
end

function plot_vorticity(sim)
	@inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
	contourf(sim.flow.σ',
			 color=palette(:BuGn), clims=(-5, 5), linewidth=0,
			 aspect_ratio=:equal, legend=false, border=:none)
end

sim = circle(32);
sim_step!(sim,10);
plot_vorticity(sim)