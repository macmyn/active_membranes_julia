using DrWatson
@quickactivate

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets

@parameters x t
@variables c(..) v(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

x_min = t_min = 0
x_max = 1
t_max = 10

c0(x,t) = 1
v0(x,t) = 1

conc_f(c) = c./(1+c) 

# eq = [Dx(c(x,t)) ~ v(x,t) - Dxx(v(x,t)),
#       Dt(c(x,t)) ~ Dxx(c(x,t)) - Dx(v(x,t)*c(x,t))]

eq = [Dx(v(x,t)) ~ 2*Dxx(v(x,t)),
      Dx(c(x,t)) ~ 2*Dxx(c(x,t))]

domains = [x ∈ Interval(x_min, x_max),
           t ∈ Interval(t_min, t_max)]

bcs = [
       c(x,0) ~ c0(x,t),
       v(x,0) ~ v0(x,t)]

@named pdesys = PDESystem(eq,bcs,domains,[x,t],[c(x,t),v(x,t)])

N = 100

order = 2

discretisation = MOLFiniteDifference([x=>N], t, approx_order=order)

println("Discretisation:")
@time prob = discretize(pdesys,discretisation)
