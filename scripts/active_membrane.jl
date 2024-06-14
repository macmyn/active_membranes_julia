using DrWatson
@quickactivate

using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Revise, Plots

@parameters x t
@variables c(..) v(..) fc(..) cv(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

x_min = t_min = 0
x_max = 1
t_max = 100

LL = 1/(2*pi)
Pe = 200

c0(x,t) = sin.(x*2*pi).^2
v0(x,t) = c0(x,t)
fc0(x,t) = c0(x,t)/(1.0.+c0(x,t))
cv0(x,t) = c0(x,t)*v0(x,t)

conc_f(c) = c./(1+c) 

# eq = [Dx(c(x,t)) ~ v(x,t) - Dxx(v(x,t)),
#       Dt(c(x,t)) ~ Dxx(c(x,t)) - Dx(v(x,t)*c(x,t))]

eq = [fc(x,t) ~ (c(x,t)/(1+c(x,t))),
      cv(x,t) ~ v(x,t) * c(x,t),
      v(x,t) ~ LL^2 * Dxx(v(x,t)) + LL*Dx(fc(x,t)),
      Dt(c(x,t)) ~ LL*Dxx(c(x,t)) - Pe*Dx(cv(x,t)),
      ]

# eq = [Dx(v(x,t)) ~ 2*Dxx(v(x,t)),
#       Dx(c(x,t)) ~ 2*Dxx(c(x,t))]

domains = [x ∈ Interval(x_min, x_max),
           t ∈ Interval(t_min, t_max)]

bcs = [
       c(x,0) ~ c0(x,t),
       v(x,0) ~ v0(x,t),  # it doesn't like it if this line isn't there...
       c(x_min,t) ~ c(x_max,t),
       Dx(c(x_min,t)) ~ Dx(c(x_max,t)),
       v(x_min,t) ~ v(x_max,t),
       Dx(v(x_min,t)) ~ Dx(v(x_max,t)),
      #  fc(x,0) ~ fc0(x,t),
      #  cv(x,0) ~ cv0(x,t)
      #  v(x_max,t) ~ 0
       ]

@named pdesys = PDESystem(eq,bcs,domains,[x,t],[c(x,t),v(x,t), fc(x,t), cv(x,t)])

N = 100

order = 2

discretisation = MOLFiniteDifference([x=>N], t, approx_order=order)

println("Discretisation:")
@time prob = discretize(pdesys,discretisation)

sol = solve(prob,Rosenbrock23())

sol_xs = sol[x]
sol_ts = sol[t]
sol_cs = sol[c(x,t)]
sol_vs = sol[v(x,t)]

anim = @animate for i in 1:length(sol_ts)
      plot(sol_xs, sol_cs[:,i],ylim=(0,1))
      plot!(sol_xs, sol_vs[:,i],ylim=(0,1))
end

gif(anim, "methodoflines.gif",fps=10)