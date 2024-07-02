using DifferentialEquations, Revise, Interpolations
using RecursiveArrayTools
include("big_eq2.jl")

# Geometric parameters %%%%%%%%%%%%%%%%%%%%%%%%%
ri = 1e-3;                     # minimum pole radius to avoid singularity 

# Geometry of the cell
R0 = 3;                        # radius of the undeformed spherical cell
vol0 = 2/3*R0^3;               # Volume of the cell

# Geometry of the constraining potential
a0 = 9/1.1;
b0 = 6/1.1;
b = 3.0;         #  The geometry of the constraint is given by
a = 4.5;         #  z^2/a^2 + r^2/b^2 = 1
sp = 10;         #  The potential of the constraint is given by
Vext = 0;        #  V = Vext*exp(sp*d), where d = z^2/a^2 + r^2/b^2 - 1 

#%%%%%%%%% Viscosity and friction %%%%%%%%%%%%%%%%%%%%%%
mu = 1;                        # viscosity mu = etas
mub = 0;                       # viscosity mub = (etab-etas)/2
Ga1 = 0;                       # normal friction coefficient
Ga2 = 0.1;                    # tangential friction coefficient

#%%%%%%%%% Motor dynamics and active stress %%%%%%%%%%%
sigeq = 1;                     # Turnover the motors is described by
ksig = 1;                      # Source term -ksig*(sigma-sigeq)
xi = 20;                       # Active stress ta = xi*sigma^2/(sigma^2+Sig0^2)
Sig0 = 1;                      # Hill-function parameter
Di = 1;                        # Diffusion constant

#%%%%%%%% Elastic properties of the surface %%%%%%%%%%%
ka = 0.1;                      # Bending rigidity (has to be > 0)
la = 1;                        # Passive membrane tension

#%%%%%%%% bvp5c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTol = 1e-3;              # relative tolerence 
Nmax = 50000;             # maximum iteration steps
dt = 0.01;                # simulation time step
T = 6;                    # total simulation time
date = "15-Jan-2024";




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
init = guessf(uc,ka,R0,ri,sigeq)
init = ArrayPartition(init...)
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
    global pert .+= ran[n]*cos.(n*uc*pi)
    global dpert .+= - ran[n]*n*pi*sin.(n*uc*pi)
end
rr = 0.2*sigeq/maximum(abs.(pert))
pert .= pert.*rr                 # This is to ensure that the maximum perturbation is no greater than 0.2*sigeq        
dpert .= dpert.*rr
init.x[9] .= init.x[9] .+ pert # Add peturbation to concentration
init.x[10] .= init.x[10] .+ dpert # Add derivative perturbation to concentration derivative


uc0 = uc
p0c = init.x[1]
r0c = init.x[3]
z0c = init.x[4]
h0c = init.x[12]


sig0c = init.x[9]

params = Dict(:press0 => press0,
              :uc0 => uc0,
              :p0c => p0c,
              :r0c => r0c,
              :z0c => z0c,
              :h0c => h0c,
              :sig0c => sig0c)

#################################################################

sols = []

bvp = BVProblem(shapef!, twobcf!, init, (0,1), params)
                                        # ^ this was uc

sol = solve(bvp, Vern8(), dt=dt)                                        
