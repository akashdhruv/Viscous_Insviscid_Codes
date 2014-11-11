function [theta] = momentum_thickness(const,s,u,du_ds)

% This function implements the momentum thickness equations derived by Katz
% and Plotking P470 and for the momentum thickness at the stagnation point
% Moran P209
% const = [cbar, rho, g, mu, nu, Re, Alpha, Delta, Vel]
% s = arc length
% u = velocity
% a = velocity gradient

rho= const(2);
mu = const(4);
nu = const(5);
% Re = const(6);

% p1 = (0.45*nu)./(u.^6);
% 
% theta(1) = sqrt((0.075*mu)/(rho*du_ds(1)));
% p1(1) = theta(1);
% 
% p2 = cumtrapz(s,u.^5);
% 
% theta = sqrt(p1.*p2);

theta = zeros(1,length(s));
theta(1) = sqrt((0.075*mu)/(rho*du_ds(1)));

in2 = cumtrapz(s,u.^5);

for i = 2:length(s)
    in1 = (0.45*nu)/(u(i)^6);
    theta(i) = sqrt(in1*in2(i));
end


% for i = 2:length(s)
%   K = 0.45/Re;
%   xm = (s(i)+s(i-1))/2;
%   dx = (s(i)-s(i-1));
%   coeff = sqrt(3/5);
% 
%   f1 = ppval(pchip(s,u),xm-coeff*dx/2); f1 = f1^5;
%   f2 = ppval(pchip(s,u),xm);            f2 = f2^5;
%   f3 = ppval(pchip(s,u),xm+coeff*dx/2); f3 = f3^5;
% 
%   dt2_du6 = K*dx/18*(5*f1+8*f2+5*f3);
%   theta(i) = theta(1)+sqrt((theta(i-1).^2*u(i-1).^6 + dt2_du6)./u(i).^6);
% end
theta = theta';

end