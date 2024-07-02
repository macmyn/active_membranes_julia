using DifferentialEquations
# using Interpolations
using Plots
using Dierckx

xs = 0.0:0.01:2π
ys = sin.(xs).^2
ll = 1/(2π) 
Pe = 10
f_c = ys./(1.0.+ys)

spl = Spline1D(xs,f_c)
df = derivative(spl, xs)
using DataInterpolations
dspl = CubicSpline(df,xs,extrapolate=true)
myf(x) = dspl(x)
# dspl = Spline1D(xs,df)
# itp = interpolate(f_c)
# spl = LinearInterpolation(xs, f_c, extrapolation_bc=Line())
# itp = interpolate((xs,),f_c,Gridded(Linear()))
# ditp = Interpolations.gradient.(itp, xs)

function simple_ode!(du, u, p, x)
    du[1] = u[2]
    du[2] = (1/ll^2) * (dspl(x) - ll*u[1])
end

u0 = [0.0, -1.0]  # Initial conditions for y and dy/dx

tspan = (0.0, 2π)

prob = ODEProblem(simple_ode!, u0, tspan)

sol = solve(prob, Rosenbrock23(), abstol=1e-10, reltol=1e-10)

y_sol = [sol(x)[1] for x in xs]

plot(xs, ys)
plot!(xs, y_sol)
plot!(xs, f_c)

# dont' know how to do periodic BCs with this. best example seems to be https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/
# but this seems unnecessarily manual.

# hence trying to use fenics