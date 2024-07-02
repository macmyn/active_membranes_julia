using ModelingToolkit
using DifferentialEquations
using MethodOfLines
using Revise
using DomainSets

@variables x t
@variables c(..) v(..)

∂t = Differential(t)
∂x = Differential(x)

# advection-diffusion constants
R = 1.0
D = 0.1
# v = 1.0
cₒ = 1.0
vₒ = 1.0
x_min = 0.0
x_max = 1.0
t_min = 0.0
t_max = 1.0
N = 200
ll = 1/(2*pi)
P = 50
xs = LinRange(0.0,1.0,N)

# c_normal(x) = exp(-(5*x)^2)
v_init(x) = 1.0
c_normal(x) = sin(((pi)*x))^2
c_init_plot(x) = sin((pi/200)*x)^2
fc(c) = c/(1+c)
# prod(x,t) = v(x,t)*c(x,t)

# @register_symbolic prod(x,t)

# advection-diffusion equation
# eqn = R * ∂t(c(x, t)) ~ D * ∂x(c(x, t))^2 - v * ∂x(c(x, t))
eqn = [∂t(c(x,t)) ~ ll*∂x(c(x,t))^2 - P*∂x(v(x,t))*c(x,t) - P*∂x(c(x,t))*v(x,t),
       (ll)*(1/(1+c(x,t)^2))*∂x(c(x,t)) ~ v(x,t) - (ll^2)*∂x(v(x,t))^2]

# initial and boundary conditions
bcs = [c(x,t_min) ~ c_normal(x),
       c(x_min,t) ~ c(x_max,t),
       v(x_min,t) ~ v(x_max,t),
       v(x,t_min) ~ v_init(x)]#,
    #    ∂x(c(x_min,t)) ~ ∂x(c(x_max,t))]

# space and time domains
dom = [x ∈ Interval(x_min, x_max), t ∈ Interval(t_min, t_max)]

# define PDE system
@named sys = PDESystem(eqn, bcs, dom, [x, t], [c(x, t),v(x,t)])

# convert the PDE into an ODE problem
prob = discretize(sys, MOLFiniteDifference([x => 200], t, advection_scheme=WENOScheme()))

# solve the problem
sol = solve(prob, Rosenbrock23(), saveat=0.01)

c = sol[c(x,t)]

using Plots

plot(xlim=(0,200),ylim=(0,10))
n = size(c, 2)
anim = @animate for t ∈ 1:n
    plot(c_init_plot,0,200)
    plot!(c[:,t],legend=false)
end
gif(anim, "diff3.gif", fps=10)