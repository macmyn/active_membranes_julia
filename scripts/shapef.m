function dy = shapef(u,y,para,ka,la,ksig,sigeq,mu,mub,Ga1,Ga2,Sig0,xi,Di,uc,r0c,z0c,p0c,h0c,sig0c,a,b,Vext,Sp,dt)

% FUNCTION SETS UP RIGHT-HAND-SIDE EQUATIONS OF ODE SYSTEM SOLVED BY bvp5c

% Read out current field values
fn0 = para(2); % Pressure
%dy = zeros(12,1);    % a column vector
p = y(1,:);    % psi
dp = y(2,:);
r = y(3,:);
z = y(4,:);
al = y(5,:);
be = y(6,:);
vu = y(7,:);
dvu = y(8,:);
sig = y(9,:);
dsig = y(10,:);
%vol = y(11);
h = y(12,:);

if Vext ~= 0
    % Constraining forces in direction r and z
    dd = z.^2./a.^2+    r.^2./b.^2-1;
    fr = -2.*r.*Sp.*Vext.*exp(Sp.*dd)./b.^2;
    fz = -2.*z.*Sp.*Vext.*exp(Sp.*dd)./a.^2;
    
    % External flow-independent forces in normal and tangential direction
    fn = fn0 +     fz.*cos(p)+    fr.*sin(p);
    fu = -fz.*sin(p)+    fr.*cos(p);
else
    % No contraining potential
    fn = fn0;
    fu = zeros(size(fn));
end

% Ga2 = 1./(1+exp(-dd));
% Ga2 = exp(dd);

% Interpolate relevant fields of previous time point solution on points u 
% requested by bvp5c solver in current step (needed to compute Euler 
% terms ~(X-X0)/dt in the bvp4c equations
r0 = spline(uc,r0c,u);
z0 = spline(uc,z0c,u);
p0 = spline(uc,p0c,u);
h0 = spline(uc,h0c,u);
sig0 = spline(uc,sig0c,u);

% Spontaneous curvature
c0 = 0;

%%% COPY HERE RIGHT-HAND SIDE EQUATIONS PRODUCED BY MATHEMATICA SCRIPT %%% 

dy(1,:)=dp;

dy(2,:)=[(-1).*dp.*h.*r.^(-1).*cos(p)+ be.*h.^2.*ka.^(-1).*r.^(-1).*cos(p)+  ...
  al.*h.^2.*ka.^(-1).*r.^(-1).*sin(p)+ h.^2.*r.^(-2).*cos(p).*sin(p) ...
  ];

dy(3,:)=h.*cos(p);

dy(4,:)=(-1).*h.*sin(p);

dy(5,:)=[(-1).*c0.*dp.*ka+ (1/2).*dp.^2.*h.^(-1).*ka+ (1/2).*c0.^2.*h.*ka+ ( ...
  -1).*fn.*h.*r.*sin(p)+ 2.*dp.*dvu.*h.^(-1).*mu.*r.*sin(p)+ 2.*dp.* ...
  dvu.*h.^(-1).*mub.*r.*sin(p)+ dp.*la.*r.*sig.^2.*(sig.^2+ Sig0.^2) ...
  .^(-1).*sin(p)+ dp.*la.*r.*Sig0.^2.*(sig.^2+ Sig0.^2).^(-1).*sin(p)+  ...
  dp.*r.*sig.^2.*(sig.^2+ Sig0.^2).^(-1).*xi.*sin(p)+ 2.*dp.*mub.*vu.* ...
  cos(p).*sin(p)+ dt.^(-1).*Ga1.*h.*r.*z.*cos(p).*sin(p)+ 2.*dp.^2.* ...
  dt.^(-1).*h.^(-1).*mu.*r.*z.*cos(p).*sin(p)+ 2.*dp.^2.*dt.^(-1).* ...
  h.^(-1).*mub.*r.*z.*cos(p).*sin(p)+ (-1).*dt.^(-1).*Ga1.*h.*r.*z0.* ...
  cos(p).*sin(p)+ (-2).*dp.^2.*dt.^(-1).*h.^(-1).*mu.*r.*z0.*cos(p).* ...
  sin(p)+ (-2).*dp.^2.*dt.^(-1).*h.^(-1).*mub.*r.*z0.*cos(p).*sin(p)+  ...
  2.*dvu.*mub.*sin(p).^2+ (-1/2).*h.*ka.*r.^(-2).*sin(p).^2+ dt.^(-1) ...
  .*Ga1.*h.*r.^2.*sin(p).^2+ 2.*dp.^2.*dt.^(-1).*h.^(-1).*mu.*r.^2.* ...
  sin(p).^2+ 2.*dp.^2.*dt.^(-1).*h.^(-1).*mub.*r.^2.*sin(p).^2+ (-1).* ...
  dt.^(-1).*Ga1.*h.*r.*r0.*sin(p).^2+ (-2).*dp.^2.*dt.^(-1).*h.^(-1) ...
  .*mu.*r.*r0.*sin(p).^2+ (-2).*dp.^2.*dt.^(-1).*h.^(-1).*mub.*r.* ...
  r0.*sin(p).^2+ h.*la.*sig.^2.*(sig.^2+ Sig0.^2).^(-1).*sin(p).^2+ h.* ...
  la.*Sig0.^2.*(sig.^2+ Sig0.^2).^(-1).*sin(p).^2+ h.*sig.^2.*(sig.^2+  ...
  Sig0.^2).^(-1).*xi.*sin(p).^2+ 2.*h.*mu.*r.^(-1).*vu.*cos(p).*sin( ...
  p).^2+ 2.*h.*mub.*r.^(-1).*vu.*cos(p).*sin(p).^2+ 4.*dp.*dt.^(-1).* ...
  mub.*z.*cos(p).*sin(p).^2+ (-4).*dp.*dt.^(-1).*mub.*z0.*cos(p).* ...
  sin(p).^2+ 4.*dp.*dt.^(-1).*mub.*r.*sin(p).^3+ (-4).*dp.*dt.^(-1).* ...
  mub.*r0.*sin(p).^3+ 2.*dt.^(-1).*h.*mu.*r.^(-1).*z.*cos(p).*sin(p) ...
  .^3+ 2.*dt.^(-1).*h.*mub.*r.^(-1).*z.*cos(p).*sin(p).^3+ (-2).*dt.^( ...
  -1).*h.*mu.*r.^(-1).*z0.*cos(p).*sin(p).^3+ (-2).*dt.^(-1).*h.* ...
  mub.*r.^(-1).*z0.*cos(p).*sin(p).^3+ 2.*dt.^(-1).*h.*mu.*sin(p).^4+  ...
  2.*dt.^(-1).*h.*mub.*sin(p).^4+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).* ...
  r0.*sin(p).^4+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*r0.*sin(p).^4];

dy(6,:)=[(-1).*fn.*h.*r.*cos(p)+ 2.*dp.*dvu.*h.^(-1).*mu.*r.*cos(p)+ 2.*dp.* ...
  dvu.*h.^(-1).*mub.*r.*cos(p)+ dp.*la.*r.*sig.^2.*(sig.^2+ Sig0.^2) ...
  .^(-1).*cos(p)+ dp.*la.*r.*Sig0.^2.*(sig.^2+ Sig0.^2).^(-1).*cos(p)+  ...
  dp.*r.*sig.^2.*(sig.^2+ Sig0.^2).^(-1).*xi.*cos(p)+ 2.*dp.*mub.*vu.* ...
  cos(p).^2+ dt.^(-1).*Ga1.*h.*r.*z.*cos(p).^2+ 2.*dp.^2.*dt.^(-1).* ...
  h.^(-1).*mu.*r.*z.*cos(p).^2+ 2.*dp.^2.*dt.^(-1).*h.^(-1).*mub.*r.* ...
  z.*cos(p).^2+ (-1).*dt.^(-1).*Ga1.*h.*r.*z0.*cos(p).^2+ (-2).* ...
  dp.^2.*dt.^(-1).*h.^(-1).*mu.*r.*z0.*cos(p).^2+ (-2).*dp.^2.*dt.^( ...
  -1).*h.^(-1).*mub.*r.*z0.*cos(p).^2+ 2.*dvu.*mub.*cos(p).*sin(p)+  ...
  dt.^(-1).*Ga1.*h.*r.^2.*cos(p).*sin(p)+ 2.*dp.^2.*dt.^(-1).*h.^(-1) ...
  .*mu.*r.^2.*cos(p).*sin(p)+ 2.*dp.^2.*dt.^(-1).*h.^(-1).*mub.* ...
  r.^2.*cos(p).*sin(p)+ (-1).*dt.^(-1).*Ga1.*h.*r.*r0.*cos(p).*sin(p) ...
  + (-2).*dp.^2.*dt.^(-1).*h.^(-1).*mu.*r.*r0.*cos(p).*sin(p)+ (-2).* ...
  dp.^2.*dt.^(-1).*h.^(-1).*mub.*r.*r0.*cos(p).*sin(p)+ h.*la.* ...
  sig.^2.*(sig.^2+ Sig0.^2).^(-1).*cos(p).*sin(p)+ h.*la.*Sig0.^2.*( ...
  sig.^2+ Sig0.^2).^(-1).*cos(p).*sin(p)+ h.*sig.^2.*(sig.^2+ Sig0.^2) ...
  .^(-1).*xi.*cos(p).*sin(p)+ 2.*h.*mu.*r.^(-1).*vu.*cos(p).^2.*sin( ...
  p)+ 2.*h.*mub.*r.^(-1).*vu.*cos(p).^2.*sin(p)+ 4.*dp.*dt.^(-1).* ...
  mub.*z.*cos(p).^2.*sin(p)+ (-4).*dp.*dt.^(-1).*mub.*z0.*cos(p).^2.* ...
  sin(p)+ 4.*dp.*dt.^(-1).*mub.*r.*cos(p).*sin(p).^2+ (-4).*dp.*dt.^( ...
  -1).*mub.*r0.*cos(p).*sin(p).^2+ 2.*dt.^(-1).*h.*mu.*r.^(-1).*z.* ...
  cos(p).^2.*sin(p).^2+ 2.*dt.^(-1).*h.*mub.*r.^(-1).*z.*cos(p).^2.* ...
  sin(p).^2+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*z0.*cos(p).^2.*sin(p) ...
  .^2+ (-2).*dt.^(-1).*h.*mub.*r.^(-1).*z0.*cos(p).^2.*sin(p).^2+ 2.* ...
  dt.^(-1).*h.*mu.*cos(p).*sin(p).^3+ 2.*dt.^(-1).*h.*mub.*cos(p).* ...
  sin(p).^3+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*r0.*cos(p).*sin(p).^3+ ( ...
  -2).*dt.^(-1).*h.*mub.*r.^(-1).*r0.*cos(p).*sin(p).^3];

dy(7,:)=dvu;

dy(8,:)=[(-1/2).*fu.*h.^2.*(mu+ mub).^(-1)+ (1/2).*Ga2.*h.^2.*(mu+ mub).^(-1) ...
  .*vu+ dsig.*h.*(mu+ mub).^(-1).*sig.^3.*(sig.^2+ Sig0.^2).^(-2).*xi+ ( ...
  -1).*dsig.*h.*(mu+ mub).^(-1).*sig.*(sig.^2+ Sig0.^2).^(-1).*xi+ (-1) ...
  .*dvu.*h.*mu.*(mu+ mub).^(-1).*r.^(-1).*cos(p)+ (-1).*dvu.*h.*mub.*( ...
  mu+ mub).^(-1).*r.^(-1).*cos(p)+ (-1).*dp.^2.*dt.^(-1).*mu.*(mu+ mub) ...
  .^(-1).*r.*cos(p)+ (-1).*dp.^2.*dt.^(-1).*mub.*(mu+ mub).^(-1).*r.* ...
  cos(p)+ dp.^2.*dt.^(-1).*mu.*(mu+ mub).^(-1).*r0.*cos(p)+ dp.^2.* ...
  dt.^(-1).*mub.*(mu+ mub).^(-1).*r0.*cos(p)+ h.^2.*mu.*(mu+ mub).^(-1) ...
  .*r.^(-2).*vu.*cos(p).^2+ h.^2.*mub.*(mu+ mub).^(-1).*r.^(-2).*vu.* ...
  cos(p).^2+ (-1).*be.*dt.^(-1).*h.^2.*ka.^(-1).*mu.*(mu+ mub).^(-1).* ...
  r.^(-1).*z.*cos(p).^2+ dp.*dt.^(-1).*h.*mub.*(mu+ mub).^(-1).*r.^( ...
  -1).*z.*cos(p).^2+ (-1).*be.*dt.^(-1).*h.^2.*ka.^(-1).*mub.*(mu+  ...
  mub).^(-1).*r.^(-1).*z.*cos(p).^2+ be.*dt.^(-1).*h.^2.*ka.^(-1).* ...
  mu.*(mu+ mub).^(-1).*r.^(-1).*z0.*cos(p).^2+ (-1).*dp.*dt.^(-1).*h.* ...
  mub.*(mu+ mub).^(-1).*r.^(-1).*z0.*cos(p).^2+ be.*dt.^(-1).*h.^2.* ...
  ka.^(-1).*mub.*(mu+ mub).^(-1).*r.^(-1).*z0.*cos(p).^2+ (-1).*dp.* ...
  dt.^(-1).*h.*mub.*(mu+ mub).^(-1).*r.^(-1).*z.*cos(2.*p)+ dp.*dt.^( ...
  -1).*h.*mub.*(mu+ mub).^(-1).*r.^(-1).*z0.*cos(2.*p)+ dp.*h.*mub.*( ...
  mu+ mub).^(-1).*r.^(-1).*vu.*sin(p)+ dp.^2.*dt.^(-1).*mu.*(mu+ mub) ...
  .^(-1).*z.*sin(p)+ dp.^2.*dt.^(-1).*mub.*(mu+ mub).^(-1).*z.*sin(p)+  ...
  (-1).*dp.^2.*dt.^(-1).*mu.*(mu+ mub).^(-1).*z0.*sin(p)+ (-1).* ...
  dp.^2.*dt.^(-1).*mub.*(mu+ mub).^(-1).*z0.*sin(p)+ dp.*dt.^(-1).*h.* ...
  mu.*(mu+ mub).^(-1).*cos(p).*sin(p)+ (-1).*be.*dt.^(-1).*h.^2.*ka.^( ...
  -1).*mu.*(mu+ mub).^(-1).*cos(p).*sin(p)+ dp.*dt.^(-1).*h.*mub.*(mu+  ...
  mub).^(-1).*cos(p).*sin(p)+ (-1).*be.*dt.^(-1).*h.^2.*ka.^(-1).* ...
  mub.*(mu+ mub).^(-1).*cos(p).*sin(p)+ (-1).*dp.*dt.^(-1).*h.*mu.*( ...
  mu+ mub).^(-1).*r.^(-1).*r0.*cos(p).*sin(p)+ be.*dt.^(-1).*h.^2.* ...
  ka.^(-1).*mu.*(mu+ mub).^(-1).*r.^(-1).*r0.*cos(p).*sin(p)+ (-1).* ...
  dp.*dt.^(-1).*h.*mub.*(mu+ mub).^(-1).*r.^(-1).*r0.*cos(p).*sin(p)+  ...
  be.*dt.^(-1).*h.^2.*ka.^(-1).*mub.*(mu+ mub).^(-1).*r.^(-1).*r0.* ...
  cos(p).*sin(p)+ (-1).*al.*dt.^(-1).*h.^2.*ka.^(-1).*mu.*(mu+ mub).^( ...
  -1).*r.^(-1).*z.*cos(p).*sin(p)+ (-1).*al.*dt.^(-1).*h.^2.*ka.^(-1) ...
  .*mub.*(mu+ mub).^(-1).*r.^(-1).*z.*cos(p).*sin(p)+ al.*dt.^(-1).* ...
  h.^2.*ka.^(-1).*mu.*(mu+ mub).^(-1).*r.^(-1).*z0.*cos(p).*sin(p)+  ...
  al.*dt.^(-1).*h.^2.*ka.^(-1).*mub.*(mu+ mub).^(-1).*r.^(-1).*z0.* ...
  cos(p).*sin(p)+ dt.^(-1).*h.^2.*mu.*(mu+ mub).^(-1).*r.^(-2).*z.* ...
  cos(p).^2.*sin(p)+ dt.^(-1).*h.^2.*mub.*(mu+ mub).^(-1).*r.^(-2).* ...
  z.*cos(p).^2.*sin(p)+ (-1).*al.*dt.^(-1).*h.^2.*ka.^(-1).*mu.*(mu+  ...
  mub).^(-1).*sin(p).^2+ (-1).*al.*dt.^(-1).*h.^2.*ka.^(-1).*mub.*( ...
  mu+ mub).^(-1).*sin(p).^2+ al.*dt.^(-1).*h.^2.*ka.^(-1).*mu.*(mu+  ...
  mub).^(-1).*r.^(-1).*r0.*sin(p).^2+ al.*dt.^(-1).*h.^2.*ka.^(-1).* ...
  mub.*(mu+ mub).^(-1).*r.^(-1).*r0.*sin(p).^2+ (-1).*dt.^(-1).*h.^2.* ...
  mu.*(mu+ mub).^(-1).*r.^(-1).*cos(p).*sin(p).^2+ (-1).*dt.^(-1).* ...
  h.^2.*mub.*(mu+ mub).^(-1).*r.^(-1).*cos(p).*sin(p).^2+ (-1/2).*dp.* ...
  dt.^(-1).*h.*mu.*(mu+ mub).^(-1).*sin(2.*p)+ (-1).*dp.*dt.^(-1).*h.* ...
  mub.*(mu+ mub).^(-1).*sin(2.*p)+ (1/2).*dp.*dt.^(-1).*h.*mu.*(mu+  ...
  mub).^(-1).*r.^(-1).*r0.*sin(2.*p)+ dp.*dt.^(-1).*h.*mub.*(mu+ mub) ...
  .^(-1).*r.^(-1).*r0.*sin(2.*p)+ (-1/2).*dt.^(-1).*h.^2.*mu.*(mu+  ...
  mub).^(-1).*r.^(-2).*z.*cos(p).*sin(2.*p)+ (-1/2).*dt.^(-1).*h.^2.* ...
  mub.*(mu+ mub).^(-1).*r.^(-2).*z.*cos(p).*sin(2.*p)+ (1/2).*dt.^(-1) ...
  .*h.^2.*mu.*(mu+ mub).^(-1).*r.^(-1).*sin(p).*sin(2.*p)+ (1/2).* ...
  dt.^(-1).*h.^2.*mub.*(mu+ mub).^(-1).*r.^(-1).*sin(p).*sin(2.*p)+  ...
  dp.*dt.^(-1).*h0.*mu.*(mu+ mub).^(-1).*sin(p+ (-1).*p0)+ dp.*dt.^(-1) ...
  .*h0.*mub.*(mu+ mub).^(-1).*sin(p+ (-1).*p0)+ dt.^(-1).*h.*h0.*mub.*( ...
  mu+ mub).^(-1).*r.^(-1).*sin(p).*sin(p+ (-1).*p0)];

dy(9,:) = dsig;

dy(10,:) = [Di.^(-1).*dvu.*h.*sig+ Di.^(-1).*dt.^(-1).*h.^2.*sig+ Di.^(-1).* ...
  h.^2.*ksig.*sig+ (-1).*Di.^(-1).*dt.^(-1).*h.^2.*sig0+ (-1).*Di.^( ...
  -1).*h.^2.*ksig.*sigeq+ Di.^(-1).*dsig.*h.*vu+ (-1).*dsig.*h.*r.^( ...
  -1).*cos(p)+ (-1).*Di.^(-1).*dsig.*dt.^(-1).*h.*r.*cos(p)+ Di.^(-1) ...
  .*dsig.*dt.^(-1).*h.*r0.*cos(p)+ Di.^(-1).*h.^2.*r.^(-1).*sig.*vu.* ...
  cos(p)+ Di.^(-1).*dp.*dt.^(-1).*h.*sig.*z.*cos(p)+ (-1).*Di.^(-1).* ...
  dp.*dt.^(-1).*h.*sig.*z0.*cos(p)+ Di.^(-1).*dp.*dt.^(-1).*h.*r.* ...
  sig.*sin(p)+ (-1).*Di.^(-1).*dp.*dt.^(-1).*h.*r0.*sig.*sin(p)+ Di.^( ...
  -1).*dsig.*dt.^(-1).*h.*z.*sin(p)+ (-1).*Di.^(-1).*dsig.*dt.^(-1).* ...
  h.*z0.*sin(p)+ Di.^(-1).*dt.^(-1).*h.^2.*r.^(-1).*sig.*z.*cos(p).* ...
  sin(p)+ (-1).*Di.^(-1).*dt.^(-1).*h.^2.*r.^(-1).*sig.*z0.*cos(p).* ...
  sin(p)+ Di.^(-1).*dt.^(-1).*h.^2.*sig.*sin(p).^2+ (-1).*Di.^(-1).* ...
  dt.^(-1).*h.^2.*r.^(-1).*r0.*sig.*sin(p).^2];

dy(11,:) = (1/2).*h.*r.^2.*sin(p);

dy(12,:) = 0;
