using DifferentialEquations, Revise, Interpolations
using StaticArrays
using BoundaryValueDiffEq
using BenchmarkTools
using Infiltrator
using RecursiveArrayTools
using ODEInterface
using FileIO,JLD2
include("big_eq2.jl")

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
# pert .= pert.*rr
pert = [   -0.1467
-0.1407
-0.1234
-0.0960
-0.0612
-0.0217
 0.0188
 0.0570
 0.0896
 0.1136
 0.1271
 0.1287
 0.1183
 0.0968
 0.0657
 0.0277
-0.0141
-0.0565
-0.0959
-0.1292
-0.1539
-0.1680
-0.1704
-0.1612
-0.1412
-0.1119
-0.0757
-0.0355
 0.0059
 0.0454
 0.0804
 0.1088
 0.1289
 0.1401
 0.1423
 0.1364
 0.1236
 0.1059
 0.0854
 0.0644
 0.0449
 0.0289
 0.0176
 0.0119
 0.0119
 0.0173
 0.0272
 0.0402
 0.0548
 0.0693
 0.0819
 0.0912
 0.0960
 0.0957
 0.0898
 0.0787
 0.0628
 0.0431
 0.0207
-0.0031
-0.0270
-0.0497
-0.0704
-0.0881
-0.1025
-0.1131
-0.1201
-0.1236
-0.1238
-0.1213
-0.1162
-0.1089
-0.0995
-0.0882
-0.0749
-0.0596
-0.0424
-0.0231
-0.0021
 0.0204
 0.0436
 0.0669
 0.0891
 0.1091
 0.1255
 0.1371
 0.1427
 0.1414
 0.1326
 0.1162
 0.0925
 0.0623
 0.0270
-0.0117
-0.0517
-0.0908
-0.1267
-0.1572
-0.1805
-0.1950
-0.2000]                 # This is to ensure that the maximum perturbation is no greater than 0.2*sigeq        
# dpert .= dpert.*rr
dpert = [         0
1.1851
2.2679
3.1548
3.7692
4.0581
3.9964
3.5891
2.8708
1.9025
0.7666
-0.4409
-1.6183
-2.6669
-3.4995
-4.0477
-4.2676
-4.1431
-3.6866
-2.9383
-1.9614
-0.8370
0.3434
1.4858
2.5016
3.3151
3.8695
4.1310
4.0913
3.7672
3.1982
2.4429
1.5724
0.6644
-0.2041
-0.9635
-1.5574
-1.9466
-2.1123
-2.0567
-1.8015
-1.3855
-0.8600
-0.2836
0.2832
0.7840
1.1710
1.4094
1.4795
1.3779
1.1173
0.7239
0.2346
-0.3069
-0.8546
-1.3644
-1.7981
-2.1264
-2.3308
-2.4044
-2.3515
-2.1864
-1.9304
-1.6094
-1.2508
-0.8801
-0.5187
-0.1822
0.1206
0.3881
0.6239
0.8362
1.0349
1.2294
1.4262
1.6272
1.8283
2.0188
2.1824
2.2979
2.3423
2.2927
2.1298
1.8407
1.4211
0.8776
0.2283
-0.4977
-1.2614
-2.0165
-2.7127
-3.3001
-3.7325
-3.9725
-3.9942
-3.7863
-3.3535
-2.7167
-1.9117
-0.9869
-0.0000]
# @infiltrate
guess_init[:,9] .= guess_init[:,9] .+ pert # Add peturbation to concentration
guess_init[:,10] .= guess_init[:,10] .+ dpert # Add derivative perturbation to concentration derivative
save_object("solve_big_eq_guess_init.jld2",guess_init)
# @infiltrate
uc0 = uc
p0c = guess_init[:,1]
r0c = guess_init[:,3]
z0c = guess_init[:,4]
h0c = guess_init[:,12]


sig0c = guess_init[:,9]

params = Dict(:ri => ri,
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
              :press0 => press0,
              :uc0 => uc0,
              :p0c => p0c,
              :r0c => r0c,
              :z0c => z0c,
              :h0c => h0c,
              :sig0c => sig0c,
              :dt => dt)

params = NamedTuple(p for p in params)


#################################################################
# @infiltrate
sols = []

bvp = BVProblem(shapef!, twobcf!, guess_init, (0,1), params)
                                        # ^ this was uc
# @infiltrate
# bvp = TwoPointBVProblem(shapef!, (bc1!,bc2!), guess_init, (0,1), params; bcresid_prototype = (zeros(9), zeros(5)))
# bcresid_prototype = (zeros(9), zeros(5)), 
println("bvp defined")
# sol = solve(bvp, Shooting(Vern8()), dt=dt)
# sol = solve(bvp, BVPM2(), dt=dt)
# sol = solve(bvp, MIRK4(), dt=dt, adaptive=false)
global dy = similar(guess_init)
shapef!(dy, guess_init, params, uc0)
println("run shapef! once")
println("solving...")
sol = solve(bvp, MIRK4(), dt=dt, adaptive=false)
println("finished")  

end


