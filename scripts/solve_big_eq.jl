using DifferentialEquations, Revise, Interpolations
using StaticArrays
using BoundaryValueDiffEq
using BenchmarkTools
using Infiltrator
using RecursiveArrayTools
using ODEInterface
using FileIO,JLD2
include("big_eq2.jl")
include("jac.jl")

# Geometric parameters %%%%%%%%%%%%%%%%%%%%%%%%%
const ri = 1e-3;                     # minimum pole radius to avoid singularity 

# Geometry of the cell
const R0 = 1;                        # radius of the undeformed spherical cell
const vol0 = 2/3*R0^3;               # Volume of the cell

# Geometry of the constraining potential
const a0 = 9/1.1;
const b0 = 6/1.1;
const b = 3.0;         #  The geometry of the constraint is given by
const a = 4.5;         #  z^2/a^2 + r^2/b^2 = 1
const sp = 10;         #  The potential of the constraint is given by
const Vext = 0;        #  V = Vext*exp(sp*d), where d = z^2/a^2 + r^2/b^2 - 1 

#%%%%%%%%% Viscosity and friction %%%%%%%%%%%%%%%%%%%%%%
const mu = 1;                        # viscosity mu = etas
const mub = 0;                       # viscosity mub = (etab-etas)/2
const Ga1 = 0;                       # normal friction coefficient
const Ga2 = 0.0;                    # tangential friction coefficient

#%%%%%%%%% Motor dynamics and active stress %%%%%%%%%%%
const sigeq = 1;                     # Turnover the motors is described by
const ksig = 0;                      # Source term -ksig*(sigma-sigeq)
const xi = 15;                       # Active stress ta = xi*sigma^2/(sigma^2+Sig0^2)
const Sig0 = 1;                      # Hill-function parameter
const Di = 2.7;                        # Diffusion constant

#%%%%%%%% Elastic properties of the surface %%%%%%%%%%%
const ka = 0.01;                      # Bending rigidity (has to be > 0)
const la = 1;                        # Passive membrane tension

#%%%%%%%% bvp5c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTol = 1e-3;              # relative tolerence 
Nmax = 50000;             # maximum iteration steps
dt = 0.01;                # simulation time step
T = 6;                    # total simulation time
date = "15-Jan-2024";

apply_pert = true

function main()


### (2) By-hand initial geometry and concentration fields
### Option 1: Spherical shape, homogeneous concentration
n_init = 101; # Number of grid points along initial circle
uc = Array(LinRange(0,1,n_init));

# Initial guess function for spherical surface
### Defined in functions file ###

# Generate bvpc compatible initial guess data
fa = xi*sigeq^2/(sigeq^2+ Sig0^2);       # Active stress
press0 = 2*(la+fa)/R0;                   # Pressure
# init = guessf(uc,ka,R0,ri,sigeq)
# init = 6
# init =[uc*pi,pi,ri.+ R0*sin.(uc*pi),R0*(cos.(uc*pi)),ka/2*ri*(0.0.-1/R0^2),0,0,0,sigeq,0,2/3*R0^3*(2 .+ cos.(uc*pi)).*sin.(uc*pi/2).^4,R0*pi]
guess_init = guessf(uc,ka,R0,ri,sigeq)
guess_init_single = guessfsingle(uc,ka,R0,ri,sigeq)
# guess_init = MArray{Tuple{101,12}}(guess_init)
# init = ArrayPartition(init...)
# My understanding is that this doesn't actually solve anything, just
# puts things in the right structure to be solved 
# sol = bvpinit(uc,yinit,[pi*R0,press0]);
#################################################################

# Apply small random perturbations to the motor concentration
ran = 0.2*rand(1,10).-0.1;
pert = zeros(size(uc))
dpert = zeros(size(uc))

# This global is currently a workaround because we're in the global scope from the REPL
for n in 1:10
    pert .+= ran[n]*cos.(n*uc*pi)
    dpert .+= - ran[n]*n*pi*sin.(n*uc*pi)
end
rr = 0.2*sigeq/maximum(abs.(pert))
pert .= pert.*rr
                 # This is to ensure that the maximum perturbation is no greater than 0.2*sigeq        
dpert .= dpert.*rr
# @infiltrate
guess_init[:,9] .= guess_init[:,9] .+ pert # Add peturbation to concentration
guess_init[:,10] .= guess_init[:,10] .+ dpert # Add derivative perturbation to concentration derivative
save_object("solve_big_eq_guess_init.jld2",guess_init)
uc0 = uc
p0c = guess_init[:,1]
r0c = guess_init[:,3]
z0c = guess_init[:,4]
h0c = guess_init[:,12]

@infiltrate

sig0c = guess_init[:,9]

params = Dict(:R0 => R0,
              :ri => ri,
              :vol0 => vol0,
              :ka => ka,
              :la => la,
              :Vext => Vext,
              :mu => mu,
              :mub => mub,
              :Ga1 => Ga1,
              :Ga2 => Ga2,
              :sigeq => sigeq,
              :ksig => ksig,
              :xi => xi,
              :Sig0 => Sig0,
              :Di => Di,
              :h => pi*R0,
              :press0 => press0,
              :uc0 => uc0,
              :p0c => p0c,
              :r0c => r0c,
              :z0c => z0c,
              :h0c => h0c,
              :sig0c => sig0c,
              :apply_pert => apply_pert,
              :dt => dt)

params = NamedTuple(p for p in params)


#################################################################
# @infiltrate
sols = []
# bvpf = BVPFunction(shapefsingle!,jacsingle)
# bvp = BVProblem(shapefsingle!, twobcfsingle!, guess_init_single, (0,1), params)
                                        # ^ this was uc
bvp = BVProblem(shapefsingle!, twobcfsingle!, initialGuess, (0,1), params)
@infiltrate
# bvp = TwoPointBVProblem(shapef!, (bc1!,bc2!), guess_init, (0,1), params; bcresid_prototype = (zeros(9), zeros(5)))
# bcresid_prototype = (zeros(9), zeros(5)), 
println("bvp defined")
# sol = solve(bvp, Shooting(Vern8()), dt=dt)
# sol = solve(bvp, BVPM2(), dt=dt)
# sol = solve(bvp, MIRK4(), dt=dt, adaptive=false)
# global dy = similar(guess_init)
# shapef!(dy, guess_init, params, uc0)
# println("run shapef! once")
println("solving...")
sol = solve(bvp, MIRK4(), dt=dt, adaptive=false)#,abstol=1e-10,reltol=1e-10)
# sol = solve(bvp, Shooting(AutoTsit5(Rosenbrock23())),dt=dt,adaptive=false)
# sol = solve(bvp, Shooting(Vern8()), dt=dt, adaptive=false)
@infiltrate
println("finished")  

return sol

end


