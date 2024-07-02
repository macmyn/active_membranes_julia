using ModelingToolkit
using MethodOfLines
using OrdinaryDiffEq
using DomainSets
using Dierckx
using DifferentialEquations

@parameters x
@variables v(..)
# Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

x_min = 0.0
x_max = 1.0
ll = 1/(2π)

xs = 0.0:0.01:2π
ys = sin.(xs).^2
f_c = ys./(1.0.+ys)

spl = Spline1D(xs,f_c)
df = derivative(spl,xs)
# dspl = Spline1D(xs,df)
using DataInterpolations
dspl = CubicSpline(df,xs)
myf(x) = dspl(x)

eq = [Dxx(v(x)) ~ (v(x) - ll * myf(x))*(1/(ll^2))]

domain = [x ∈ Interval(x_min,x_max)]

bcs = [v(x_min) ~ v(x_max)]

tspan = (x_min, x_max)

@named sys = PDESystem(eq,bcs,domain,[x],[v(x)])

N = 32
order = 2
discretisation = MOLFiniteDifference([x=>N])
prob = discretize(sys,discretisation)

sol = solve(prob, Tsit5(), saveat=0.2)
sol(xs)