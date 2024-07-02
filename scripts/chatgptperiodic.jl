using DifferentialEquations
using Plots

# Constants
const D = 1.0       # Diffusion coefficient
const L = 1.0       # Length of the domain
const Nx = 100      # Number of spatial points
const dx = L / (Nx-1)   # Spatial step size
const x = range(0, L, length=Nx)  # Spatial grid

# Initial condition: Concentration profile at t=0
function initial_condition(x)
    μ = L / 2
    σ = 0.1
    return exp(-(x .- μ).^2 / (2σ^2))
end

# Diffusion equation with periodic boundary conditions
function diffusion_eq!(du, u, p, t)
    D = p[1]
    Nx = length(u)
    for i in 1:Nx
        du[i] = D * (u[mod1(i+1, Nx)] - 2u[i] + u[mod1(i-1, Nx)]) / dx^2
    end
end

# Define the periodic boundary condition callback function
function periodic_boundary_condition!(integrator)
    u = integrator.u
    Nx = length(u)
    u[1] = u[end-1]
    u[end] = u[2]
end

# Define a condition that always returns true to apply the callback at every step
condition(u, t, integrator) = true

# Create the DiscreteCallback for periodic boundary conditions
cb = DiscreteCallback(condition, periodic_boundary_condition!)

# Initial condition and time span
u0 = initial_condition.(x)
tspan = (0.0, 1.0)

# Define the problem
prob = ODEProblem(diffusion_eq!, u0, tspan, [D])

# Solve the problem
sol = solve(prob, Tsit5(), callback=cb)

# Plot the concentration profile over time
anim = @animate for i in 1:length(sol.t)
    plot(x, sol.u[i], ylim=(0, 1), label="t=$(sol.t[i])", xlabel="x", ylabel="Concentration")
end

gif(anim, "diffusion_periodic.gif", fps=10)
