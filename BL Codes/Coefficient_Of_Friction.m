function [dstarj, thetaj, cfricj, shapej, TLTSu, TLTSl] = Coefficient_Of_Friction(const, xj, zj, sj, uej, aej, SP)

% Boundary Layer Theory Setup 
% Written by Christopher J. Blower
% Date: 2013-11-27
% Reminder: The velocity profile has to be modeled in absolute values,
% consequently allowing the Re to be positive across the upper and lower
% surface accordingly. 
% Seperate Vectors in Upper and Lower Surfaces about Stagnation Point (SP). 
% Each surface is to be defined from the stagnation point. Leading Edge to
% Trailing Edge.

%% Constants
ref   = const(12)+16;
n     = const(8)+1;                                                               % # of Points Distributed Along Each Surface
Re    = const(6);                                                          % Free-stream Reynolds Number
V     = const(7);                                                          % Free-stream Velocity
cbar  = const(1);                                                          % Chord Length
rho   = const(2);                                                          % Air Density
mu = const(4);                                                             % Dynamic Viscosity

%% Airfoil Parameters
s_lsm= linspace(sj(1),SP(3),n);                                            % Lower Surface - Sectioned Arc Length LSTE to LSLE
s_ls = flipud((SP(3)*ones(1,n))'-(s_lsm)');                                % Lower Surface [Modified] - Arc Length from LE to TE
x_ls = flipud((ppval(pchip(sj,xj),s_lsm))');                               % Lower Surface - X-Coordinate
z_ls = flipud((ppval(pchip(sj,zj),s_lsm))');                               % Lower Surface - Z-Coordinate
ue_ls= abs(flipud((ppval(pchip(sj,uej),s_lsm))'));                         % Lower Surface - Velocity Profile
ae_ls= abs(flipud((ppval(pchip(sj,aej),s_lsm))'));                         % Lower Surface - Velocity Gradient

s_usm= linspace(SP(3),sj(end),n);                                          % Upper Surface - Sectioned Arc Length USLE to USTE
s_us = (s_usm-(SP(3)*ones(1,n)))';                                         % Upper Surface [Modified] - Arc Length from LE to TE
x_us = (ppval(pchip(sj,xj),s_usm))';                                       % Upper Surface - X-Coordinate
z_us = (ppval(pchip(sj,zj),s_usm))';                                       % Upper Surface - Z-Coordinate
ue_us= ((ppval(pchip(sj,uej),s_usm))');                                    % Upper Surface - Velocity Profile
ae_us= (ppval(pchip(sj,aej),s_usm))';                                      % Upper Surface - Velocity Gradient

% figure; 
% plot(x_us,ae_us)

ue_us(1) = 0;                                                              % Define Stagnation Point Velocity [Upper Surface]
ue_ls(1) = 0;                                                              % Define Stagnation Point Velocity [Lower Surface]  

USP = [x_us, z_us, s_us, ue_us, ae_us];                                    % Upper Surface Profile Data - X- , Z-Coordinates, Arc length, Velocity, Velocity Gradient
LSP = [x_ls, z_ls, s_ls, ue_ls, ae_ls];                                    % Lower Surface Profile Data - X- , Z-Coordinates, Arc length, Velocity, Velocity Gradient
 
%% Boundary Layer Solver
[x_us, x_ls, s_blus, s_blls, del_1_us, del_1_ls, del_2_us, del_2_ls, cfric_us, cfric_ls,H_us, H_ls, TLTSu, TLTSl] = boundary_layer_solver(const, LSP, USP);

TLTSu;
TLTSl;

% Note: The boundary layer solver implements 500 nodes over the airfoil
% surface, thereby allowing the position of the transition point to be
% accurate than that of an airfoil with only 30 per surface. 
ind = (length(xj)+1)/2;
%% Redistribute to Panel End Points
dstar_u = ppval(pchip(s_blus,del_1_us),s_us);                              % Displacement Thickness - Upper Surface
dstar_l = ppval(pchip(s_blls,del_1_ls),s_ls);                              % Displacement Thickness - Lower Surface

theta_u = ppval(pchip(x_us,del_2_us),flipud(xj(1:ind)));                    % Momentum Thickness - Upper Surface
theta_l = ppval(pchip(x_ls,del_2_ls),flipud(xj(1:ind)));                    % Momentum Thickness - Lower Surface

% figure;
% plot(s_us,theta_u,'-r',s_ls,theta_l,'-b')

cfric_u = ppval(pchip(s_blus,cfric_us),s_us);                              % Coefficient of Friction - Upper Surface
cfric_l = ppval(pchip(s_blls,cfric_ls),s_ls);                              % Coefficient of Friction - Lower Surface

shape_u = ppval(pchip(s_blus,H_us),s_us);                                  % Shape Factor - Upper Surface
shape_l = ppval(pchip(s_blls,H_ls),s_ls);                                  % Shape Factor - Lower Surface

%% Define Boundary Parameter Vectors Clockwise [LSTE-LE-USTE] 
dstarj = [flipud(dstar_l(2:end)); dstar_u];                                % Displacement Thickness at End Points LSTE-LE-USTE
thetaj = [flipud(theta_l(2:end)); theta_u];                                % Momentum Thickness at End Points LSTE-LE-USTE
cfricj = [flipud(cfric_l(2:end)); cfric_u];                                % Coefficient of Friciton at End Points LSTE-LE-USTE
shapej = [flipud(shape_l(2:end)); shape_u];                                % Shape Factor at End Points LSTE-LE-USTE
       
% Cd = sy(res_BLu(1),res_BLu(2),res_BLu(3),res_BLl(1),res_BLl(2),res_BLl(3));                                                          % Coefficient of Drag
end