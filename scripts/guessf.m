function g = guessf(u,ka,R0,ri,sigeq)
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
% vol =y(11);
% h = y(12);
% ka = 1;
%ka = 1;
g = [u*pi;...
     pi;...
     ri+ R0*sin(u*pi);...
     R0*(cos(u*pi));...
     ka/2*ri*(0-1/R0^2);...
     0;...
     0;...
     0;...
     sigeq;...
     0;...
     2/3*R0^3*(2+ cos(u*pi))*sin(u*pi/2)^4;...
     R0*pi];