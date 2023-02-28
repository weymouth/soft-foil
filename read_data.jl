using Plots, CSV 

d(x) = x<0.3 ? 0 : (x-0.3)^2

plot(;xaxis="x/L",yaxis="y/L")
for pressure in ["0","0.2","0.5","0.7"]
    file = CSV.File("data/Base Foil Section "*pressure*" bar.csv")
    x = filter(!ismissing,file["x centre line"])
    y = filter(!ismissing,file[" y centre line"])
    L = maximum(x)
    plot!(x/L,y/L,label=pressure*" bar")
    x = 0:0.01:1
    A = maximum(y/L)/0.7^2
    plot!(x, A*d.(x), label="",ls=:dash)
end
savefig("centerline with fit.png")

file = CSV.File("data/Base Foil Section 0 bar.csv")
x = filter(!ismissing,file["x centre line"])
y = filter(!ismissing,file[" y centre line"])
L = maximum(x)
A = maximum(y/L)/0.7^2
x = file["x points"]/L
y = file["y points"]/L
plot(x,abs.(y-A*d.(x)))
savefig("section shape.png")
