using DifferentialEquations
using Interpolations

## Test getting v from c0

c0 = sin(x/pi)
x_min = 0.0
x_max = 1.0
nx = 100
xs = LinRange(x_min,x_max,nx)
c_init = c0.(xs)


# function v_from_c(du, u, p, t)
#     dc = only.(Interpolations.gradient.(interpolate((xs,),p[:c],Gridded(Linear())),xs))
#     v = u[1]
#     dv = u[2]
#     du[1] = dv
#     du[2] = u[1] - (1./(1+p[:c]).^2) * dc
# end

# prob = ODEProblem(v_from_c, c_init, (0.0,1.0))

function second_deriv(du,u,p,t)
    v = u[1]
    dv = u[2]
    du[1] = dv
    du[2] = p[:c]
end

prob = ODEProblem(second_deriv, (0.0,1.0))