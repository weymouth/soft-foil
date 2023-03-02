using Plots,Roots

P0,P1,P2,P3 = [0,0],[0,0.1],[0.5,0.12],[1.,0.]
A,B,C,D = P0,-3P0+3P1,3P0-6P1+3P2,-P0+3P1-3P2+P3
curve(t) = A+B*t+C*t^2+D*t^3
dcurve(t) = B+2C*t+3D*t^2

candidates(X) = union(0,1,find_zeros(t -> (X-curve(t))'*dcurve(t),0,1,naive=true,no_pts=3,xatol=0.01))
function distance(X,t)
    V = X-curve(t)
    s = sign([V[2],-V[1]]'*dcurve(t))
    return s*√sum(abs2,V)
end
sdf(X,time) = argmin(abs,distance(X,t) for t in candidates(X))

x = range(-0.5, 1.5, length=100)
y = range(0, 0.5, length=50)
f(x,y) = sdf([x,y],0)
z = @. f(x',y);
contour(x,y,z,aspect_ratio=:equal, legend=false, border=:none)
xy = curve.(0:0.01:1)
plot!(first.(xy),last.(xy),ls=:dash)
savefig("distance function.png")
