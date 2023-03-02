using Plots, CSV

plot(;xaxis="x/L",yaxis="y/L")
for pressure in ["0","0.2","0.5","0.7"]
    file = CSV.File("data/Base Foil Section "*pressure*" bar.csv")
    x = filter(!ismissing,file["x centre line"])
    y = filter(!ismissing,file[" y centre line"])
    L = maximum(x)
    plot!(x/L,y/L,label=pressure*" bar")
end
x = 0:0.01:1
d(x) = x<0.3 ? 0 : (x-0.3)^2/0.7^2
plot!(x, 0d.(x), label="d/L=0",ls=:dash)
plot!(x, 0.06d.(x), label="d/L=0.06",ls=:dash)
plot!(x, 0.12d.(x), label="d/L=0.12",ls=:dash)
savefig("centerline with fit.png")

file = CSV.File("data/Base Foil Section 0 bar.csv")
x = filter(!ismissing,file["x centre line"])
y = filter(!ismissing,file[" y centre line"])
L = maximum(x)
Amp = maximum(y/L)
x = file["x points"]/L
y = file["y points"]/L
plot(x,abs.(y-Amp*d.(x)),label="centered section")

P0,P1,P2,P3 = [0,0],[0,0.1],[0.5,0.12],[1.,0.]
scatter!(first.([P0,P1,P2,P3]),last.([P0,P1,P2,P3]),label="control points")
A,B,C,D = P0,-3P0+3P1,3P0-6P1+3P2,-P0+3P1-3P2+P3
curve(t) = A+B*t+C*t^2+D*t^3
xy = curve.(0:0.01:1)
plot!(first.(xy),last.(xy),label="cubic",ls=:dash)
savefig("section shape.png")
