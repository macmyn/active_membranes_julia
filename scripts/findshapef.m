clear all;

%%%%%%%%% Geometric parameters %%%%%%%%%%%%%%%%%%%%%%%%%
ri = 1e-3;                     % minimum pole radius to avoid singularity 

% Geometry of the cell
R0 = 3;                        % radius of the undeformed spherical cell
vol0 = 2/3*R0^3;               % Volume of the cell

% Geometry of the constraining potential
a0 = 9/1.1;
b0 = 6/1.1;
b = 3.0;         %  The geometry of the constraint is given by
a = 4.5;         %  z^2/a^2 + r^2/b^2 = 1
sp = 10;         %  The potential of the constraint is given by
Vext = 0;        %  V = Vext*exp(sp*d), where d = z^2/a^2 + r^2/b^2 - 1 

%%%%%%%%%% Viscosity and friction %%%%%%%%%%%%%%%%%%%%%%
mu = 1;                        % viscosity mu = etas
mub = 0;                       % viscosity mub = (etab-etas)/2
Ga1 = 0;                       % normal friction coefficient
Ga2 = 0.1;                    % tangential friction coefficient

%%%%%%%%%% Motor dynamics and active stress %%%%%%%%%%%
sigeq = 1;                     % Turnover the motors is described by
ksig = 1;                      % Source term -ksig*(sigma-sigeq)
xi = 20;                       % Active stress ta = xi*sigma^2/(sigma^2+Sig0^2)
Sig0 = 1;                      % Hill-function parameter
Di = 1;                        % Diffusion constant

%%%%%%%%% Elastic properties of the surface %%%%%%%%%%%
ka = 0.1;                      % Bending rigidity (has to be > 0)
la = 1;                        % Passive membrane tension

%%%%%%%%% bvp5c parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTol = 1e-3;              % relative tolerence 
Nmax = 50000;             % maximum iteration steps
dt = 0.01;                % simulation time step
T = 6;                    % total simulation time
date = '15-Jan-2024';

rng('shuffle');

%%%% Load the steady state solution as the initial condition for the dynamics
% DIR = '/Users/amietke/Uni/Oxford/Research/ActiveSurfaces/IsotropicCode/Unconstrained/dynamic/data/';
% filename = [DIR,'Shape_ConstraintPotential_SteadyState',...
%    '_la_',num2str(la),...
%    '_R0_',num2str(R0),...
%    '_mu_',num2str(mu),...
%    '_mub_',num2str(mub),...
%    '_Ga1_',num2str(Ga1),...
%    '_Ga2_',num2str(Ga2),...
%    '_xi_',num2str(0),...
%    '_Di_',num2str(Di),...
%    '_Sig0_',num2str(Sig0),...
%    '_sigeq_',num2str(sigeq),...
%    '_ksig_',num2str(ksig),...
%    '_a0_',num2str(a0),...
%    '_b0_',num2str(b0),...
%    '_Vext_',num2str(Vext),...
%    '_Sp_',num2str(sp),...
%    '_RTol_',num2str(RTol),...
%    '_date_',date1,'.mat'];

%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITION SELECTION %%%%%%%%%%%%%%%%%%%%%%%%
%%% (1) Load from another solution, e.g. constrained shape initial conditions
% filename = ['/Users/amietke/Uni/Oxford/Research/ActiveSurfaces/IsotropicCode/',... 
%     'ConstrainingPotential/dynamic/InitData/Shape_ConstraintPotential_dynamics_',...
%     'la_1_R0_3_mu_1_mub_0_Ga1_0_Ga2_1_xi_0_Di_1_Sig0_1_sigeq_1_ksig_1_a_4.5_b_3_',...
%     'Vext_0.5_Sp_10_RTol_0.001_date_07-Jan-2024.mat'];
% 
% load(filename);
% sol = sol_all{1,2};

%%% (2) By-Hand initial geometry and concentration fields
%%% Option 1: Spherical shape, homogeneous concentration
n_init = 101; % Number of grid points along initial circle
uc = linspace(0,1,n_init);

% Initial guess function for spherical surface
yinit = @(u) guessf(u,ka,R0,ri,sigeq);

% Generate bvpc compatible initial guess data
fa = xi*sigeq^2/(sigeq^2+ Sig0^2);       % Active stress
press0 = 2*(la+fa)/R0;                   % Pressure
sol = bvpinit(uc,yinit,[pi*R0,press0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply small perturbations to the motor concentration
ran = 0.2*rand(1,10)-0.1;
pert = 0;
dpert = 0;

% Grid coordinates of the solution
uc = sol.x;

%%%%%% Add concentration perturbation %%%%%%
for n = 1:10
    pert = pert+ran(n)*cos(n*uc*pi);
    dpert = dpert - ran(n)*n*pi*sin(n*uc*pi);
end
rr = 0.2*sigeq/max(abs(pert));
pert = pert*rr;                 % This is to ensure that the maximum perturbation is no greater than 0.2*sigeq        
dpert = dpert*rr;
sol.y(9,:) = sol.y(9,:) + pert; % Add peturbation to concentration
sol.y(10,:) = sol.y(10,:) + dpert; % Add derivative perturbation to concentration derivative
   
% Save the data in the following file
DIR = pwd;
filename = [DIR,'/SimData/Shape_dynamics',...
   '_ka_',num2str(ka),...
   '_la_',num2str(la),...
   '_R0_',num2str(R0),...
   '_mu_',num2str(mu),...
   '_mub_',num2str(mub),...
   '_Ga1_',num2str(Ga1),...
   '_Ga2_',num2str(Ga2),...
   '_xi_',num2str(xi),...
   '_Di_',num2str(Di),...
   '_Sig0_',num2str(Sig0),...
   '_sigeq_',num2str(sigeq),...
   '_ksig_',num2str(ksig),...
   '_a_',num2str(a),...
   '_b_',num2str(b),...
   '_Vext_',num2str(Vext),...
   '_Sp_',num2str(sp),...
   '_RTol_',num2str(RTol),...
   '_dt_',num2str(dt),...
   '_T_',num2str(T),...
   '_date_',date,'.mat'];

sol_all = cell(100,4); % Empty array to store intermediate solutions
Tr = linspace(0,T,100); % Time points to save data
j = 0; % Step counter
for t = 0:dt:T
    % General content of the solution vector
    % [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, 12] = ...
    % [Psi, Psi', r, z, alpha, beta, vu, vu', sig, sig', Volume, h]

    % Lagrange multupliers for global constraints 
    % para = [h0 = Total arc length ,Pressure]
    para = sol.parameters;    
   
    % Current reference grid point coordinates and geometry
    uc0 = sol.x;
    p0c = sol.y(1,:); % Tangent angle psi
    r0c = sol.y(3,:); % r(u)
    z0c = sol.y(4,:); % z(u)    
    h0c = sol.y(12,:); % Coordinate scaling factor h0

    % Current reference concentration
    sig0c = sol.y(9,:);
    
    % Current reference solution
    sol0 = sol;
    
    % Input functions for bvp5c solver
    yeq = @(u,y,para) shapef(u,y,para,ka,la,ksig,sigeq,mu,mub,Ga1,Ga2,Sig0,xi,Di,uc,r0c,z0c,p0c,h0c,sig0c,a0,b0,Vext,sp,dt);
    ybc = @(ya,yb,para) twobcf(ya,yb,para,ka,sp,ri,vol0,Ga2);
    jach = @(u,y,para) jac(u,y,para,ka,la,ksig,sigeq,mu,mub,Ga1,Ga2,Sig0,xi,Di,uc,r0c,z0c,p0c,h0c,sig0c,a0,b0,Vext,sp,dt);
    opts = bvpset('RelTol',RTol,'AbsTol',1e-8,'NMax',Nmax,'FJacobian',jach,'Vectorized','on');
    
    % Determine bvp5c solution for current shape and concentration
    sol = bvp5c(yeq,ybc,sol,opts);
    uc = sol.x;
    yc = deval(sol,uc);
    pc = yc(1,:);
    rc = yc(3,:);
    zc = yc(4,:);        
    
    % Get reference r(u) and z(u) from previous reference solution at new grid points
    if strcmp(sol0.solver,'bvpinit')
        % If previous data came from a by-hand guess/initial conditions
        r0c = spline(uc0,r0c,uc);
        z0c = spline(uc0,z0c,uc);
    elseif (strcmp(sol0.solver,'bvp4c'))||(strcmp(sol0.solver,'bvp5c'))
        % If previous data comes from a bvpc solution                               
        y0 = deval(sol0,uc);
        r0c = y0(3,:);
        z0c = y0(4,:);
    end
    
    % dr/dt and dz/dt
    dr = (rc - r0c)/dt; 
    dz = (zc - z0c)/dt;

    % Tangent angle and local normal velocity
    psi = pc;
    vn = cos(psi).*dz + sin(psi).*dr;

    if min(abs(Tr-t)) < 0.5*dt
        % Save solution at specific time steps
        j = j+1;
        sol_all{j,1} = t;
        sol_all{j,2} = sol;
        sol_all{j,3} = vn;
        save(filename,'sol_all') 
        
        % Show maximal normal velocity
        disp(max(abs(vn)));

        % Show current time step
        disp(t);        
    end
    
    if sol.stats.maxerr > RTol
        % Abort if required error tolerance 
        % is not satisfied anymore
        break;
    end
    
end



