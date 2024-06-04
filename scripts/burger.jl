using DomainSets
using MethodOfLines
using ModelingToolkit
using OrdinaryDiffEq

@parameters x t
@variables u(..), u_sq(..)
Dx = Differential(x)
Dt = Differential(t)

x_min = 0.0
x_max = 1.0
t_min = 0.0
t_max = 1.0

eq =  [u_sq(t,x) ~ u(t,x)^2, 
       Dt(u(t, x)) ~ -0.5 * Dx(u_sq(t,x))]

bcs = [u(0, x) ~ (1 + sin(2*π*x)),
       u(t, x_min) ~ u(t, x_max),
       u_sq(t, x_min) ~ u_sq(t, x_max)]

domains = [t ∈ Interval(t_min, t_max),
    x ∈ Interval(x_min, x_max)]

dx = 0.001

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x), u_sq(t, x)])

disc = MOLFiniteDifference([x => dx], t)

prob = discretize(pdesys, disc)