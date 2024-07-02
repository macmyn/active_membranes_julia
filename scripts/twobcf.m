function res = twobcf(ya,yb,para,ka,sp,ri,vol0,Ga2)

% psi = y(1);
% dpsi = y(2);
% r = y(3);
% z = y(4);
% alpha = y(5);
% beta = y(6);
% vu = y(7);
% dvu = y(8);
% sigma = y(9);
% dsigma = y(10);
% vol = y(11);
% h = y(12);
% ka = 1;

c0 = 0;
h = para(1);
alpha_i = ka/2*ri*(c0^2-ya(2)^2/h^2);

if Ga2 == 0 % No tangential friction
    beta_z_bc0 = ya(6) + yb(6); % beta(0) + beta(1) = 0
    beta_z_bc1 = ya(4) + yb(4); % z(0) + z(1) = 0
else
    beta_z_bc0 = ya(6); % beta(1) = 0
    beta_z_bc1 = yb(6); % beta(1) = 0
end

res = [ ya(1)-0;...   % psi(0) = 0;
        ya(3)-ri;...      % r(0) = ri;
        beta_z_bc0;...       
        ya(5)-alpha_i;... % alpha(0)=ka/2*(c0^2-dpsi^2/h^2)...
        beta_z_bc1;...      
        ya(7)-0;...       % vu(0) = 0;
        ya(10)-0;...      % sigma'(0)=0;
        ya(11)-0;...      % vol(0) = 0;
        ya(12)-h;...      % h(0) = h0;
        yb(1)-pi;...      % psi(1) = pi;
        yb(3)-ri;...      % r(1) = ri;
        yb(7)-0;...       % vu(1) = 0;
        yb(10)-0;...      % sigma'(1) = 0;
        yb(11)-vol0];     % vol(1) = vol0;