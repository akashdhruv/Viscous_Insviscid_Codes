function [rho,c,mu,nu,Cp,Cv,lambda]=air_properties(temperature,pressure)
% Computes air density and speed of sound as a function of temperature and pressure
% [rho,c,[mu,nu,Cp,Cv,lambda]]=FM_cstes(temperature, [pressure]);
% inputs :
% Symbol            name                  units
% temperature    temperature             Celsius
% optional :
% pressure     atmospheric pressure        Pa
%
% returns :
% Symbol            name                  units 
% Rho            air density              kg.m-3
% c             speed of sound            m.s-1
% optional :
% mu           dynamic viscosity        kg.m-1.s-1
% nu          kinematic viscosity         m2.s-1
% Cp     specific heat constant pressure
% Cv      specific heat constant volume
% lambda   coefficient of thermal Conductivity    W.m-1.K-1
% input parameter pressure is optional default is atmospheric pressure

	% Cp/Cv
gamma = 1.400;
	% Characteristic temperature
T_e = 78.6;
	% Molecular diameter
sigma = 3.711;
	% Temp in K
T = 273.15+temperature;
	% atmospheric pressure by default
if nargin==1,
	pressure= 101300;
end;
	% Molecular mass
Mw=28.94;
	% Boltzmann gas constant
Ru = 8315;

rho = pressure*Mw/(Ru*T);
c = sqrt(gamma*pressure/rho);

if nargout>=3,
Omega_v = 1.147*(T/T_e)^(-0.145) + (T/T_e+0.5)^(-2);
mu = 26.69e-7*sqrt(Mw*T)/(sigma^2*Omega_v);
end;

if nargout>=4,
nu = mu/rho;
end;

if nargout>=5,
Cv = (gamma-1)*Ru/Mw;
end;
if nargout>=6,
Cp = Cv*gamma;
end;
if nargout>=7,
%k0 = 0.0241;
%lambda = k0*(T/273)^0.81;
lambda = 5/3*mu*Cv;
end;


%EOF
