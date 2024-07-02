function simple_ode!(du, u, p, x)
    du .= cos(x)
end

u0 = [1.0] # initial condition
tspan = (0.0, 2π) # time span for the solution

prob = ODEProblem(simple_ode!, u0, tspan)
sol = solve(prob, abstol=1e-10, reltol=1e-10)

xs = 0.0:0.01:2π
ys = [sol(x)[1] for x in xs]

plot(xs, ys)