function [x_us, x_ls, s_us, s_ls DSTAR_u, DSTAR_l, THETA_u, THETA_l, CFRIC_u, CFRIC_l, SHAPE_u, SHAPE_l, TTLSu, TTLSl] = boundary_layer_solver(const, LSP, USP)

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
nu    = const(5);                                                          % Kinematic Viscosity of Air 
V     = const(7);                                                          % Freestream Velocity [m/s]
n     = const(8)+1;                                                      % # of Points Distributed Along Each Surface
Re    = const(6);                                                          % Freestream Reynolds Number
TSH = 2.7;                                                                 % Turbulent Shape Factor Separation
%% Airfoil Parameters

s_ls = LSP(:,3);                                                           % Lower Surface - Arc Length from LE to TE
x_ls = LSP(:,1);                                                           % Lower Surface - X-Coordinate
z_ls = LSP(:,2);                                                           % Lower Surface - Z-Coordinate
ue_ls= LSP(:,4);                                                           % Lower Surface - Velocity Profile
ae_ls= LSP(:,5);                                                           % Lower Surface - Velocity Gradient

s_us = USP(:,3);                                                           % Upper Surface - Arc Length from LE to TE
x_us = USP(:,1);                                                           % Upper Surface - X-Coordinate
z_us = USP(:,2);                                                           % Upper Surface - Z-Coordinate
ue_us= USP(:,4);                                                           % Upper Surface - Velocity Profile
ae_us= USP(:,5);                                                           % Upper Surface - Velocity Gradient

%% Boundary Layer Parameters
Re_xu = abs((ue_us.*x_us)/nu);                                             % Reynolds Number - X-Position Dependant [Upper Surface]
Re_xl = abs((ue_ls.*x_ls)/nu);                                             % Reynolds Number - X-Position Dependant [Lower Surface]

Re_su = abs((ue_us.*s_us)/nu);                                             % Reynolds Number - Arc Length Dependant [Upper Surface]
Re_sl = abs((ue_ls.*s_ls)/nu);                                             % Reynolds Number - Arc Length Dependant [Lower Surface]

[del_2_Lu] = momentum_thickness(const,s_us,ue_us,ae_us);                   % Momentum Thickness for Steady, 2D, Incompressible LAMINAR Flow [Katz P470]
[del_2_Ll] = momentum_thickness(const,s_ls,ue_ls,ae_ls);                   % Momentum Thickness for Steady, 2D, Incompressible LAMINAR Flow [Katz P470]

Lambda_u = ((del_2_Lu.^2)/nu).*ae_us;                                      % Dimensionless Pressure-Gradient Parameter [Upper Surface]
Lambda_l = ((del_2_Ll.^2)/nu).*ae_ls;                                      % Dimensionless Pressure-Gradient Parameter [Lower Surface]

% Lambda_u = Re*(del_2_Lu.^2).*ae_us;                                        % Dimensionless Pressure-Gradient Parameter [Upper Surface]
% Lambda_l = Re*(del_2_Ll.^2).*ae_ls;                                        % Dimensionless Pressure-Gradient Parameter [Lower Surface]

%% Thwaites' Parameters        
[~, H_Lu] = Cebeci_n_Bradshaw(Lambda_u);                                   % Cebeci & Bradshaw Parameters - Laminar [Upper Surface]
[~, H_Ll] = Cebeci_n_Bradshaw(Lambda_l);                                   % Cebeci & Bradshaw Parameters - Laminar [Lower Surface]

[del_1_Lu] = displacement_thickness(H_Lu,del_2_Lu);                        % Displacement Thickness for Steady, 2D, Incompressible LAMINAR Flow [Upper Surface]
[del_1_Ll] = displacement_thickness(H_Ll,del_2_Ll);                        % Displacement Thickness for Steady, 2D, Incompressible LAMINAR Flow [Lower Surface]

Re_d2_u = abs(((ue_us).*del_2_Lu)/nu);                                     % Reynolds Number Momentum Thickness Dependent [Upper Surface]
Re_d2_l = abs(((ue_ls).*del_2_Ll)/nu);                                     % Reynolds Number Momentum Thickness Dependent [Lower Surface]

[cf_Lu] = laminar_cf(del_1_Lu, del_2_Lu, ue_us, s_us);                     % Laminar Coefficient of Friction - Upper Surface
[cf_Ll] = laminar_cf(del_1_Ll, del_2_Ll, ue_ls, s_ls);                     % Laminar Coefficient of Friction - Upper Surface

% figure;
% plot(x_us,H_Lu,'-b',...
%      x_ls,H_Ll,'-r')
% axis([0 1 1 3.5])
%% Transition Point - Michels Criterion [Upper & Lower Surfaces]
[x_tpu,z_tpu,s_tpu] = Michel_Transition(Re_su,Re_d2_u,x_us,z_us,s_us);     % Michels Transition Point - Upper Surface
[x_tpl,z_tpl,s_tpl] = Michel_Transition(Re_sl,Re_d2_l,x_ls,z_ls,s_ls);     % Michels Transition Point - Lower Surface

%         figure;
%         plot(x_us,H_u,'-b',...
%              x_ls,H_l,'-r',...
%              [x_tpu x_tpu], [0 4],'--b',...
%              [x_tpl x_tpl], [0 4],'--r')

LAS_usx = linspace(x_us(1),x_tpu,n);                                       % Laminar Attached Surface Flow - x-coord - Upper Surface     
LAS_usz = ppval(pchip(x_us,z_us),LAS_usx);                                 % Laminar Attached Surface Flow - z-coord - Upper Surface
LAS_uss = ppval(pchip(x_us,s_us),LAS_usx);                                 % Laminar Attached Surface Flow - s-coord - Upper Surface

LAS_lsx = linspace(x_ls(1),x_tpl,n);                                       % Laminar Attached Surface Flow - x-coord - Lower Surface
LAS_lsz = ppval(pchip(x_ls,z_ls),LAS_lsx);                                 % Laminar Attached Surface Flow - z-coord - Lower Surface
LAS_lss = ppval(pchip(x_ls,s_ls),LAS_lsx);                                 % Laminar Attached Surface Flow - s-coord - Lower Surface

%         figure;
%         plot(LAS_usx,LAS_usz,'-b',...
%              LAS_lsx,LAS_lsz,'-r')
%         axis([0 1 -0.5 0.5])
%         axis equal 
% 
%         figure;
%         plot(x_us,del_1_Lu,'-r',...
%              x_ls,del_1_Ll,'-b')
 
TP_iu = find(ge(x_us,x_tpu),1);
TP_il = find(ge(x_ls,x_tpl),1);

if gt(x_tpu,x_us(end))
%     disp('Laminar Flow Only')
    
    THETA_u = del_2_Lu;
    THETA_l = del_2_Ll;
    DSTAR_u = del_1_Lu;
    DSTAR_l = del_1_Ll;
    CFRIC_u = cf_Lu;
    CFRIC_l = cf_Ll;
    SHAPE_u = H_Ll;
    SHAPE_l = H_Ll;
    
else
%     disp('Transition Incurred')
    %% Turbulent Modeling
    [cf_Tu,H1_tu,H_Tu,del_2_Tu,del_1_Tu] = heads_turbulence(TP_iu,x_us,s_us,ue_us,Re_d2_u,del_1_Lu,del_2_Lu,H_Lu,cf_Lu,Re);  
    [cf_Tl,H1_tl,H_Tl,del_2_Tl,del_1_Tl] = heads_turbulence(TP_il,x_ls,s_ls,ue_ls,Re_d2_l,del_1_Ll,del_2_Ll,H_Ll,cf_Ll,Re);

    %% Join Laminar & Turbulent Vectors
    THETA_u = [del_2_Lu(1:TP_iu-2);
              (del_2_Lu(TP_iu-1)+del_2_Tu(1))/2;
              (del_2_Lu(TP_iu)+del_2_Tu(2))/2;
               del_2_Tu(3:end)'];

    THETA_l = [del_2_Ll(1:TP_il-2);
              (del_2_Ll(TP_il-1)+del_2_Tl(1))/2;
              (del_2_Ll(TP_il)+del_2_Tl(2))/2;
               del_2_Tl(3:end)'];

    %         figure;
    %         plot(x_us, THETA_u,'-r',...
    %              x_ls, THETA_l,'-b')         

    DSTAR_u = [del_1_Lu(1:TP_iu-2);
              (del_1_Lu(TP_iu-1)+del_1_Tu(1))/2;
              (del_1_Lu(TP_iu)+del_1_Tu(2))/2;
               del_1_Tu(3:end)'];

    DSTAR_l = [del_1_Ll(1:TP_il-2);
              (del_1_Ll(TP_il-1)+del_1_Tl(1))/2;
              (del_1_Ll(TP_il)+del_1_Tl(2))/2;
               del_1_Tl(3:end)'];

    %         figure;
    %         plot(x_us, DSTAR_u,'-r',...
    %              x_ls, DSTAR_l,'-b')

    CFRIC_u = [cf_Lu(1:TP_iu-2)';
              (cf_Lu(TP_iu-1)+cf_Tu(1))/2;
              (cf_Lu(TP_iu)+cf_Tu(2))/2;
               cf_Tu(3:end)'];

    CFRIC_l = [cf_Ll(1:TP_il-2)';
              (cf_Ll(TP_il-1)+cf_Tl(1))/2;
              (cf_Ll(TP_il)+cf_Tl(2))/2;
               cf_Tl(3:end)'];

    %          figure;
    %          plot(x_us, CFRIC_u,'-r',...
    %               x_ls, CFRIC_l,'-b')

    SHAPE_u = [H_Lu(1:TP_iu-2)';
              (H_Lu(TP_iu-1)+H_Tu(1))/2;
              (H_Lu(TP_iu)+H_Tu(2))/2;
               H_Tu(3:end)'];

    SHAPE_l = [H_Ll(1:TP_il-2)';
              (H_Ll(TP_il-1)+H_Tl(1))/2;
              (H_Ll(TP_il)+H_Tl(2))/2;
               H_Tl(3:end)'];

    % for i = 1:length(SHAPE_l)
    %     if eq(SHAPE_l(1),0)
    %         SHAPE_l(1) = 0;
    %     elseif eq(SHAPE_l(i),0)
    %         SHAPE_l(i) = SHAPE_l(i-1);
    %     end
    % end    
    COMPP_u = [Lambda_u, DSTAR_u, THETA_u, SHAPE_u];

    % Smoothing Function
    [Lambda_u, DSTAR_u, THETA_u, SHAPE_u, CFRIC_u] = bl_smoothing(x_us,Lambda_u, DSTAR_u, THETA_u, SHAPE_u, CFRIC_u);
    [Lambda_l, DSTAR_l, THETA_l, SHAPE_l, CFRIC_l] = bl_smoothing(x_ls,Lambda_l, DSTAR_l, THETA_l, SHAPE_l, CFRIC_l);

end

%              figure;
%              plot(x_us, SHAPE_u,'-b',...
%                   x_ls, SHAPE_l,'-r')
          

%% Find Laminar Separation
lam_sep_u = find(le(Lambda_u,-0.09),1);                                    % Laminar Separation Node - Upper Surface
lam_sep_l = find(le(Lambda_l,-0.09),1);                                    % Laminar Separation Node - Lower Surface

if isempty(lam_sep_u)
    x_lsu = x_us(end);
else 
    x_lsu = x_us(lam_sep_u(1));
%     if ge(x_lsu,x_tpu)
%         x_lsu = x_us(end);
%     end
end

z_lsu = ppval(pchip(x_us,z_us),x_lsu);
s_lsu = ppval(pchip(x_us,s_us),x_lsu);

if isempty(lam_sep_l)
    x_lsl = x_ls(end);
else 
    x_lsl = x_ls(lam_sep_l(1));
%     if ge(x_lsl,x_tpl)
%         x_lsl = x_ls(end);
%     end
end

z_lsl = ppval(pchip(x_ls,z_ls),x_lsl);
s_lsl = ppval(pchip(x_ls,s_ls),x_lsl);

%% Find Turbulent Separation
x_TU = find(ge(x_us,x_tpu));
x_usT= x_us(x_TU);
SHAPE_TU = SHAPE_u(x_TU);                                                  % Shape Factor within Turbulent Region
turb_sep_u = find(ge(SHAPE_TU,TSH));                                       % Turbulent Separation Node - Upper Surface

if isempty(turb_sep_u)
    x_tsu = x_us(end);
else 
    x_tsu = x_usT(turb_sep_u(1));                                          % Turbulent Separation X-Coordinate
end

x_tsu;
z_tsu = ppval(pchip(x_us,z_us),x_tsu);
s_tsu = ppval(pchip(x_us,s_us),x_tsu);

x_TL = find(ge(x_ls,x_tpl));
x_lsT= x_ls(x_TL);
SHAPE_TL = SHAPE_l(x_TL);
turb_sep_l = find(ge(SHAPE_TL,TSH));                                       % Shape Factor of 3.0 for Turbulent Separation...

if isempty(turb_sep_l)
    x_tsl = x_ls(end);
else 
    x_tsl = x_lsT(turb_sep_l(1));                                           % Turbulent Separation X-Coordinate
end

z_tsl = ppval(pchip(x_ls,z_ls),x_tsl);
s_tsl = ppval(pchip(x_ls,s_ls),x_tsl);

%% Compile Reference Points

TTLSu = [x_tpu, z_tpu, s_tpu;                                              % Transition Point Coordinates
         x_lsu, z_lsu, s_lsu;                                              % Laminar Separation Point Coordinates         
         x_tsu, z_tsu, s_tsu];                                             % Turbulent Separation Point Coordinates

TTLSl = [x_tpl, z_tpl, s_tpl;                                              % Transition Point Coordinates
         x_lsl, z_lsl, s_lsl;                                              % Laminar Separation Point Coordinates         
         x_tsl, z_tsl, s_tsl];                                             % Turbulent Separation Point Coordinates

%% Boundary Layer Thickness - Full Airfoil
SHAPE_u(1) = 0;
SHAPE_l(1) = 0;
CFRIC_u(1) = 0;
CFRIC_l(1) = 0;
% [del_99_LTT_u] = Boundary_Layer_Profile(const, n, s_us, ue_us, del_99_FP_UL, Re_su, s_tpu_D); 
% [del_99_LTT_l] = Boundary_Layer_Profile(const, n, s_ls, ue_ls, del_99_FP_LL, Re_sl, s_tpl_D);

end