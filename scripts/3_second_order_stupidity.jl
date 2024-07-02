using DifferentialEquations
using Interpolations
using Plots


xs = 0.0:0.01:2π
ys = cos.(xs)

itp = LinearInterpolation(xs, ys, extrapolation_bc=Line())

function simple_ode!(du, u, p, x)
    du[1] = u[2]
    du[2] = itp(x)    
end

u0 = [1.0, 0.0]

tspan = (0.0, 2π)

prob = ODEProblem(simple_ode!, u0, tspan)

sol = solve(prob, abstol=1e-10, reltol=1e-10)

y_sol = [sol(x)[1] for x in xs]

plot(xs, ys)
plot!(xs, y_sol)