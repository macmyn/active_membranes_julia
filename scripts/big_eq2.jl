using Revise

function shapef!(dy, y, params, u)

# FUNCTION SETS UP RIGHT-HAND-SIDE EQUATIONS OF ODE SYSTEM SOLVED BY bvp5c
# @infiltrate
# y = y[1]  # for some reason y (ArrayPartiion) is being passed as Vecor(ArrayPartition) i.e. [y]

# Read out current field values
# fn0 = para(2); # Pressure
#dy = zeros(12,1);    # a column vector
# @infiltrate
p = y[:,1];    # psi
dp = y[:,2];
r = y[:,3];
z = y[:,4];
al = y[:,5];
be = y[:,6];
vu = y[:,7];
dvu = y[:,8];
sig = y[:,9];
dsig = y[:,10];
#vol = y(11);
h = y[:,12];

fn0 = params[:press0]

if Vext != 0
    # Constraining forces in direction r and z
    dd = z.^2.0./a.^2 +    r.^2.0./b.^2-1;
    fr = -2.0.*r.*Sp.*Vext.*exp(Sp.*dd)./b.^2;
    fz = -2.0.*z.*Sp.*Vext.*exp(Sp.*dd)./a.^2;
    
    # External flow-independent forces in normal and tangential direction
    fn = fn0 +     fz.*cos.(p)+    fr.*sin.(p);
    fu = -fz.*sin.(p)+    fr.*cos.(p);
else
    # No contraining potential
    fn = fn0;
    fu = zeros(size(fn));
end

# Ga2 = 1./(1.+exp(-dd));
# Ga2 = exp(dd);

uc0 = params[:uc0]
r0c = params[:r0c]
p0c = params[:p0c]
z0c = params[:z0c]
h0c = params[:h0c]
sig0c = params[:sig0c]
dt = params[:dt]

# Interpolate relevant fields of previous time point solution on points u 
# requested by bvp5c solver in current step (needed to compute Euler 
# terms ~(X-X0)/dt in the bvp4c equations
r0itp = interpolate((uc0,), r0c, Gridded(Linear()))
r0 = r0itp.(u)  # Interpolate at the points in u;
z0itp = interpolate((uc0,), z0c, Gridded(Linear()))
z0 = z0itp.(u)  # Interpolate at the points in u;
p0itp = interpolate((uc0,), p0c, Gridded(Linear()))
p0 = p0itp.(u)  # Interpolate at the points in u;
h0itp = interpolate((uc0,), h0c, Gridded(Linear()))
h0 = h0itp.(u)  # Interpolate at the points in u;
sig0itp = interpolate((uc0,), sig0c, Gridded(Linear()))
sig0 = sig0itp.(u)  # Interpolate at the points in u;

# Spontaneous curvature
c0 = 0;

### COPY HERE RIGHT-HAND SIDE EQUATIONS PRODUCED BY MATHEMATICA SCRIPT ### 
# println(dy)
# println("DP", dp)
# @infiltrate

dy[:,1].=dp;

dy[:,2].=(-1).*dp.*h.*r.^(-1).*cos.(p).+ be.*h.^2.0.*ka.^(-1).*r.^(-1).*cos.(p).+al.*h.^2.0.*ka.^(-1).*r.^(-1).*sin.(p).+ h.^2.0.*r.^(-2).*cos.(p).*sin.(p) 
  

dy[:,3].=h.*cos.(p);

dy[:,4].=(-1).*h.*sin.(p);

dy[:,5].=(-1).*c0.*dp.*ka.+ (1/2).*dp.^2.0.*h.^(-1).*ka.+ (1/2).*c0.^2.0.*h.*ka.+ (-1).*fn.*h.*r.*sin.(p).+ 2.0.*dp.*dvu.*h.^(-1).*mu.*r.*sin.(p).+ 2.0.*dp.*dvu.*h.^(-1).*mub.*r.*sin.(p).+ dp.*la.*r.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*sin.(p).+ dp.*la.*r.*Sig0.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*sin.(p).+ dp.*r.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*xi.*sin.(p).+ 2.0.*dp.*mub.*vu.*cos.(p).*sin.(p).+ dt.^(-1).*Ga1.*h.*r.*z.*cos.(p).*sin.(p).+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*z.*cos.(p).*sin.(p).+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*z.*cos.(p).*sin.(p).+ (-1).*dt.^(-1).*Ga1.*h.*r.*z0.*cos.(p).*sin.(p).+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*z0.*cos.(p).*sin.(p).+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*z0.*cos.(p).*sin.(p).+ 2.0.*dvu.*mub.*sin.(p).^2.0.+ (-1/2).*h.*ka.*r.^(-2).*sin.(p).^2.0.+ dt.^(-1).*Ga1.*h.*r.^2.0.*sin.(p).^2.0.+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.^2.0.*sin.(p).^2.0.+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.^2.0.*sin.(p).^2.0.+ (-1).*dt.^(-1).*Ga1.*h.*r.*r0.*sin.(p).^2.0.+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*r0.*sin.(p).^2.0.+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*r0.*sin.(p).^2.0.+ h.*la.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*sin.(p).^2.0.+ h.*la.*Sig0.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*sin.(p).^2.0.+ h.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*xi.*sin.(p).^2.0.+ 2.0.*h.*mu.*r.^(-1).*vu.*cos.(p).*sin.(p).^2.0.+ 2.0.*h.*mub.*r.^(-1).*vu.*cos.(p).*sin.(p).^2.0.+ 4.0*dp.*dt.^(-1).*mub.*z.*cos.(p).*sin.(p).^2.0.+ (-4).*dp.*dt.^(-1).*mub.*z0.*cos.(p).*sin.(p).^2.0.+ 4.0*dp.*dt.^(-1).*mub.*r.*sin.(p).^3.0.+ (-4).*dp.*dt.^(-1).*mub.*r0.*sin.(p).^3.0.+ 2.0.*dt.^(-1).*h.*mu.*r.^(-1).*z.*cos.(p).*sin.(p).^3.0.+ 2.0.*dt.^(-1).*h.*mub.*r.^(-1).*z.*cos.(p).*sin.(p).^3.0.+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*z0.*cos.(p).*sin.(p).^3.0.+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*z0.*cos.(p).*sin.(p).^3.0.+ 2.0.*dt.^(-1).*h.*mu.*sin.(p).^4.0.+ 2.0.*dt.^(-1).*h.*mub.*sin.(p).^4.0.+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*r0.*sin.(p).^4.0.+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*r0.*sin.(p).^4;

dy[:,6].=(-1).*fn.*h.*r.*cos.(p).+ 2.0.*dp.*dvu.*h.^(-1).*mu.*r.*cos.(p).+ 2.0.*dp.*dvu.*h.^(-1).*mub.*r.*cos.(p).+ dp.*la.*r.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*cos.(p).+ dp.*la.*r.*Sig0.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*cos.(p).+ dp.*r.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*xi.*cos.(p).+ 2.0.*dp.*mub.*vu.*cos.(p).^2.0.+ dt.^(-1).*Ga1.*h.*r.*z.*cos.(p).^2.0.+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*z.*cos.(p).^2.0.+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*z.*cos.(p).^2.0.+ (-1).*dt.^(-1).*Ga1.*h.*r.*z0.*cos.(p).^2.0.+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*z0.*cos.(p).^2.0.+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*z0.*cos.(p).^2.0.+ 2.0.*dvu.*mub.*cos.(p).*sin.(p).+ dt.^(-1).*Ga1.*h.*r.^2.0.*cos.(p).*sin.(p).+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.^2.0.*cos.(p).*sin.(p).+ 2.0.*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.^2.0.*cos.(p).*sin.(p).+ (-1).*dt.^(-1).*Ga1.*h.*r.*r0.*cos.(p).*sin.(p).+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mu.*r.*r0.*cos.(p).*sin.(p).+ (-2).*dp.^2.0.*dt.^(-1).*h.^(-1).*mub.*r.*r0.*cos.(p).*sin.(p).+ h.*la.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*cos.(p).*sin.(p).+ h.*la.*Sig0.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*cos.(p).*sin.(p).+ h.*sig.^2.0.*(sig.^2.0.+ Sig0.^2).^(-1).*xi.*cos.(p).*sin.(p).+ 2.0.*h.*mu.*r.^(-1).*vu.*cos.(p).^2.0.*sin.(p).+ 2.0.*h.*mub.*r.^(-1).*vu.*cos.(p).^2.0.*sin.(p).+ 4.0*dp.*dt.^(-1).*mub.*z.*cos.(p).^2.0.*sin.(p).+ (-4).*dp.*dt.^(-1).*mub.*z0.*cos.(p).^2.0.*sin.(p).+ 4.0*dp.*dt.^(-1).*mub.*r.*cos.(p).*sin.(p).^2.0.+ (-4).*dp.*dt.^(-1).*mub.*r0.*cos.(p).*sin.(p).^2.0.+ 2.0.*dt.^(-1).*h.*mu.*r.^(-1).*z.*cos.(p).^2.0.*sin.(p).^2.0.+ 2.0.*dt.^(-1).*h.*mub.*r.^(-1).*z.*cos.(p).^2.0.*sin.(p).^2.0.+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*z0.*cos.(p).^2.0.*sin.(p).^2.0.+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*z0.*cos.(p).^2.0.*sin.(p).^2.0.+ 2.0.*dt.^(-1).*h.*mu.*cos.(p).*sin.(p).^3.0.+ 2.0.*dt.^(-1).*h.*mub.*cos.(p).*sin.(p).^3.0.+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*r0.*cos.(p).*sin.(p).^3.0.+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*r0.*cos.(p).*sin.(p).^3;

dy[:,7].=dvu;

dy[:,8].=(-1/2).*fu.*h.^2.0.*(mu.+ mub).^(-1).+ (1/2).*Ga2.*h.^2.0.*(mu.+ mub).^(-1).*vu.+ dsig.*h.*(mu.+ mub).^(-1).*sig.^3.0.*(sig.^2.0.+ Sig0.^2).^(-2).*xi.+ (-1).*dsig.*h.*(mu.+ mub).^(-1).*sig.*(sig.^2.0.+ Sig0.^2).^(-1).*xi.+ (-1) .*dvu.*h.*mu.*(mu.+ mub).^(-1).*r.^(-1).*cos.(p).+ (-1).*dvu.*h.*mub.*( mu.+ mub).^(-1).*r.^(-1).*cos.(p).+ (-1).*dp.^2.0.*dt.^(-1).*mu.*(mu.+ mub) .^(-1).*r.*cos.(p).+ (-1).*dp.^2.0.*dt.^(-1).*mub.*(mu.+ mub).^(-1).*r.* cos.(p).+ dp.^2.0.*dt.^(-1).*mu.*(mu.+ mub).^(-1).*r0.*cos.(p).+ dp.^2.0.* dt.^(-1).*mub.*(mu.+ mub).^(-1).*r0.*cos.(p).+ h.^2.0.*mu.*(mu.+ mub).^(-1).*r.^(-2).*vu.*cos.(p).^2.0.+ h.^2.0.*mub.*(mu.+ mub).^(-1).*r.^(-2).*vu.* cos.(p).^2.0.+ (-1).*be.*dt.^(-1).*h.^2.0.*ka.^(-1).*mu.*(mu.+ mub).^(-1).* r.^(-1).*z.*cos.(p).^2.0.+ dp.*dt.^(-1).*h.*mub.*(mu.+ mub).^(-1).*r.^( -1).*z.*cos.(p).^2.0.+ (-1).*be.*dt.^(-1).*h.^2.0.*ka.^(-1).*mub.*(mu.+  mub).^(-1).*r.^(-1).*z.*cos.(p).^2.0.+ be.*dt.^(-1).*h.^2.0.*ka.^(-1).* mu.*(mu.+ mub).^(-1).*r.^(-1).*z0.*cos.(p).^2.0.+ (-1).*dp.*dt.^(-1).*h.* mub.*(mu.+ mub).^(-1).*r.^(-1).*z0.*cos.(p).^2.0.+ be.*dt.^(-1).*h.^2.0.* ka.^(-1).*mub.*(mu.+ mub).^(-1).*r.^(-1).*z0.*cos.(p).^2.0.+ (-1).*dp.* dt.^(-1).*h.*mub.*(mu.+ mub).^(-1).*r.^(-1).*z.*cos.(2.0.*p).+ dp.*dt.^( -1).*h.*mub.*(mu.+ mub).^(-1).*r.^(-1).*z0.*cos.(2.0.*p).+ dp.*h.*mub.*( mu.+ mub).^(-1).*r.^(-1).*vu.*sin.(p).+ dp.^2.0.*dt.^(-1).*mu.*(mu.+ mub) .^(-1).*z.*sin.(p).+ dp.^2.0.*dt.^(-1).*mub.*(mu.+ mub).^(-1).*z.*sin.(p).+  (-1).*dp.^2.0.*dt.^(-1).*mu.*(mu.+ mub).^(-1).*z0.*sin.(p).+ (-1).* dp.^2.0.*dt.^(-1).*mub.*(mu.+ mub).^(-1).*z0.*sin.(p).+ dp.*dt.^(-1).*h.* mu.*(mu.+ mub).^(-1).*cos.(p).*sin.(p).+ (-1).*be.*dt.^(-1).*h.^2.0.*ka.^( -1).*mu.*(mu.+ mub).^(-1).*cos.(p).*sin.(p).+ dp.*dt.^(-1).*h.*mub.*(mu.+  mub).^(-1).*cos.(p).*sin.(p).+ (-1).*be.*dt.^(-1).*h.^2.0.*ka.^(-1).* mub.*(mu.+ mub).^(-1).*cos.(p).*sin.(p).+ (-1).*dp.*dt.^(-1).*h.*mu.*( mu.+ mub).^(-1).*r.^(-1).*r0.*cos.(p).*sin.(p).+ be.*dt.^(-1).*h.^2.0.* ka.^(-1).*mu.*(mu.+ mub).^(-1).*r.^(-1).*r0.*cos.(p).*sin.(p).+ (-1).* dp.*dt.^(-1).*h.*mub.*(mu.+ mub).^(-1).*r.^(-1).*r0.*cos.(p).*sin.(p).+  be.*dt.^(-1).*h.^2.0.*ka.^(-1).*mub.*(mu.+ mub).^(-1).*r.^(-1).*r0.* cos.(p).*sin.(p).+ (-1).*al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mu.*(mu.+ mub).^( -1).*r.^(-1).*z.*cos.(p).*sin.(p).+ (-1).*al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mub.*(mu.+ mub).^(-1).*r.^(-1).*z.*cos.(p).*sin.(p).+ al.*dt.^(-1).* h.^2.0.*ka.^(-1).*mu.*(mu.+ mub).^(-1).*r.^(-1).*z0.*cos.(p).*sin.(p).+  al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mub.*(mu.+ mub).^(-1).*r.^(-1).*z0.* cos.(p).*sin.(p).+ dt.^(-1).*h.^2.0.*mu.*(mu.+ mub).^(-1).*r.^(-2).*z.* cos.(p).^2.0.*sin.(p).+ dt.^(-1).*h.^2.0.*mub.*(mu.+ mub).^(-1).*r.^(-2).* z.*cos.(p).^2.0.*sin.(p).+ (-1).*al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mu.*(mu.+  mub).^(-1).*sin.(p).^2.0.+ (-1).*al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mub.*( mu.+ mub).^(-1).*sin.(p).^2.0.+ al.*dt.^(-1).*h.^2.0.*ka.^(-1).*mu.*(mu.+  mub).^(-1).*r.^(-1).*r0.*sin.(p).^2.0.+ al.*dt.^(-1).*h.^2.0.*ka.^(-1).* mub.*(mu.+ mub).^(-1).*r.^(-1).*r0.*sin.(p).^2.0.+ (-1).*dt.^(-1).*h.^2.0.* mu.*(mu.+ mub).^(-1).*r.^(-1).*cos.(p).*sin.(p).^2.0.+ (-1).*dt.^(-1).* h.^2.0.*mub.*(mu.+ mub).^(-1).*r.^(-1).*cos.(p).*sin.(p).^2.0.+ (-1/2).*dp.* dt.^(-1).*h.*mu.*(mu.+ mub).^(-1).*sin.(2.0.*p).+ (-1).*dp.*dt.^(-1).*h.* mub.*(mu.+ mub).^(-1).*sin.(2.0.*p).+ (1/2).*dp.*dt.^(-1).*h.*mu.*(mu.+  mub).^(-1).*r.^(-1).*r0.*sin.(2.0.*p).+ dp.*dt.^(-1).*h.*mub.*(mu.+ mub) .^(-1).*r.^(-1).*r0.*sin.(2.0.*p).+ (-1/2).*dt.^(-1).*h.^2.0.*mu.*(mu.+  mub).^(-1).*r.^(-2).*z.*cos.(p).*sin.(2.0.*p).+ (-1/2).*dt.^(-1).*h.^2.0.* mub.*(mu.+ mub).^(-1).*r.^(-2).*z.*cos.(p).*sin.(2.0.*p).+ (1/2).*dt.^(-1).*h.^2.0.*mu.*(mu.+ mub).^(-1).*r.^(-1).*sin.(p).*sin.(2.0.*p).+ (1/2).* dt.^(-1).*h.^2.0.*mub.*(mu.+ mub).^(-1).*r.^(-1).*sin.(p).*sin.(2.0.*p).+  dp.*dt.^(-1).*h0.*mu.*(mu.+ mub).^(-1).*sin.(p.+ (-1).*p0).+ dp.*dt.^(-1).*h0.*mub.*(mu.+ mub).^(-1).*sin.(p.+ (-1).*p0).+ dt.^(-1).*h.*h0.*mub.*( mu.+ mub).^(-1).*r.^(-1).*sin.(p).*sin.(p.+ (-1).*p0);

dy[:,9] .= dsig;

dy[:,10] .= Di.^(-1).*dvu.*h.*sig.+ Di.^(-1).*dt.^(-1).*h.^2.0.*sig.+ Di.^(-1).* h.^2.0.*ksig.*sig.+ (-1).*Di.^(-1).*dt.^(-1).*h.^2.0.*sig0.+ (-1).*Di.^( -1).*h.^2.0.*ksig.*sigeq.+ Di.^(-1).*dsig.*h.*vu.+ (-1).*dsig.*h.*r.^( -1).*cos.(p).+ (-1).*Di.^(-1).*dsig.*dt.^(-1).*h.*r.*cos.(p).+ Di.^(-1) .*dsig.*dt.^(-1).*h.*r0.*cos.(p).+ Di.^(-1).*h.^2.0.*r.^(-1).*sig.*vu.* cos.(p).+ Di.^(-1).*dp.*dt.^(-1).*h.*sig.*z.*cos.(p).+ (-1).*Di.^(-1).* dp.*dt.^(-1).*h.*sig.*z0.*cos.(p).+ Di.^(-1).*dp.*dt.^(-1).*h.*r.* sig.*sin.(p).+ (-1).*Di.^(-1).*dp.*dt.^(-1).*h.*r0.*sig.*sin.(p).+ Di.^( -1).*dsig.*dt.^(-1).*h.*z.*sin.(p).+ (-1).*Di.^(-1).*dsig.*dt.^(-1).* h.*z0.*sin.(p).+ Di.^(-1).*dt.^(-1).*h.^2.0.*r.^(-1).*sig.*z.*cos.(p).* sin.(p).+ (-1).*Di.^(-1).*dt.^(-1).*h.^2.0.*r.^(-1).*sig.*z0.*cos.(p).* sin.(p).+ Di.^(-1).*dt.^(-1).*h.^2.0.*sig.*sin.(p).^2.0.+ (-1).*Di.^(-1).* dt.^(-1).*h.^2.0.*r.^(-1).*r0.*sig.*sin.(p).^2;

dy[:,11].= (1/2).*h.*r.^2.0.*sin.(p);

dy[:,12].= zeros(size(p));

end

function guessf(u,ka,R0,ri,sigeq)
  # psi = y(1);
  # dpsi = y(2);
  # r = y(3);
  # z = y(4);
  # alpha = y(5);
  # beta = y(6);
  # vu = y(7);
  # dvu = y(8);
  # sigma = y(9);
  # dsigma = y(10);
  # vol =y(11);
  # h = y(12);
  # ka = 1;
  #ka = 1;
  g = [u*pi,
       ones(size(u))*pi,
       ri.+ R0*sin.(u*pi),
       R0*(cos.(u*pi)),
       ones(size(u))*ka/2*ri*(0.0.-1/R0^2),
       ones(size(u))*0.0,
       ones(size(u))*0.0,
       ones(size(u))*0.0,
       ones(size(u))*sigeq,
       ones(size(u))*0.0,
       2/3*R0^3*(2 .+ cos.(u*pi)).*sin.(u*pi/2).^4,
       ones(size(u))*R0*pi]
  return g
end

# function twobcf!(ya,yb,para,ka,sp,ri,vol0,Ga2)
# function twobcf!(res, y, p, x)
function twobcf!(res, y, p, x)
  # psi = y(1);
  # dpsi = y(2);
  # r = y(3);
  # z = y(4);
  # alpha = y(5);
  # beta = y(6);
  # vu = y(7);
  # dvu = y(8);
  # sigma = y(9);
  # dsigma = y(10);
  # vol = y(11);
  # h = y(12);
  # ka = 1;
  c0 = 0;
#   h = para(1);
  # h = y[:,12]
  h = p[:h0c]  # althought shouldn't h be the line above...?
  y = y[1]  # for some reason, y is being passed as [y]??
  # @infiltrate

  @assert allequal(h)
  # @infiltrate
  alpha_i = ka/2*ri*(c0^2-y[:,2][begin]^2/h[1].^2);
  
  if Ga2 == 0 # No tangential friction
    #   beta_z_bc0 = ya(6) .+ yb(6); # beta(0) .+ beta(1) = 0
    #   beta_z_bc1 = ya(4) .+ yb(4); # z(0) .+ z(1) = 0
    beta_z_bc0 = y[:,6][begin] .+ y[:,6][end]
    beta_z_bc1 = y[:,4][begin] .+ y[:,4][end]
  else
    #   beta_z_bc0 = ya(6); # beta(1) = 0
    #   beta_z_bc1 = yb(6); # beta(1) = 0
    beta_z_bc0 = y[:,6][begin]
    beta_z_bc1 = y[:,6][end]
  end

  # @infiltrate

  res = [ y[:,1][begin]-0;   # psi(0) = 0;
          y[:,3][begin]-ri;      # r(0) = ri;
          beta_z_bc0;       
          y[:,5][begin]-alpha_i; # alpha(0)=ka/2*(c0^2-dpsi^2/h^2)...
          beta_z_bc1;    # 
          y[:,7][begin]-0;       # vu(0) = 0;
          y[:,10][begin]-0;      # sigma'(0)=0;
          y[:,11][begin]-0;      # vol(0) = 0;
          y[:,12][begin]-h[1];      # h(0) = h0;
          y[:,1][end]-pi;      # psi(1) = pi;
          y[:,3][end]-ri;      # r(1) = ri;
          y[:,7][end]-0;       # vu(1) = 0;
          y[:,10][end]-0;      # sigma'(1) = 0;
          y[:,11][end]-vol0];     # vol(1) = vol0;
    # res[1] = ya[1]-0;   # psi(0) = 0;
    # res[2] = ya[3]-ri;      # r(0) = ri;
    # res[3] = beta_z_bc0;       
    # res[4] = ya[5]-alpha_i; # alpha(0)=ka/2*(c0^2-dpsi^2/h^2)...
    # res[5] = beta_z_bc1;    # 
    # res[6] = ya[7]-0;       # vu(0) = 0;
    # res[7] = ya[10]-0;      # sigma'(0)=0;
    # res[8] = ya[11]-0;      # vol(0) = 0;
    # res[9] = ya[12]-h;      # h(0) = h0;
    # res[10] = yb[1]-pi;      # psi(1) = pi;
    # res[11] = yb[3]-ri;      # r(1) = ri;
    # res[12] = yb[7]-0;       # vu(1) = 0;
    # res[13] = yb[10]-0;      # sigma'(1) = 0;
    # res[14] = yb[11]-vol0;     # vol(1) = vol0;


end