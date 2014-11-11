function [cpj] = CP_extrap(xi,zi,si,cpi,xj,zj,sj,SP)
% Note: The Cp currently cannot be modelled against the arc length because
% it is currently referenced against the xi components, which is the
% midpoints of each panel rather than the end points. Consequently, the
% results have to be extrapolated to allow the end points to also be
% included. 
% Because the arc length has been removed from both the beginning and the
% end of the airfoil profile beacuse of it measuring the midpoints the
% distance cut from both the beginning and the end of the airfoil has to be
% calculated.
% Reference Fig [1.2] in PhD Book 002 P 008

si_us0 = 0;                                                                % Upper Surface Arc Point 00
si_us1 = si_us0+sqrt((xi(end-3)-xi(end-4))^2+(zi(end-3)-zi(end-4))^2);     % Upper Surface Arc Point 01
si_us2 = si_us1+sqrt((xi(end-2)-xi(end-3))^2+(zi(end-2)-zi(end-3))^2);     % Upper Surface Arc Point 02
si_us3 = si_us2+sqrt((xi(end-1)-xi(end-2))^2+(zi(end-1)-zi(end-2))^2);     % Upper Surface Arc Point 03
si_us4 = si_us3+sqrt((xi(end)-xi(end-1))^2+(zi(end)-zi(end-1))^2);         % Upper Surface Arc Point 04
si_us5 = si_us4+sqrt((xj(end)-xi(end))^2+(zj(end)-zi(end))^2);             % Upper Surface Arc Point 05

cp_us0 = cpi(end-4);                                                       % Upper Surface Cp Point 00                      
cp_us1 = cpi(end-3);                                                       % Upper Surface Cp Point 01
cp_us2 = cpi(end-2);                                                       % Upper Surface Cp Point 02
cp_us3 = cpi(end-1);                                                       % Upper Surface Cp Point 03
cp_us4 = cpi(end);                                                         % Upper Surface Cp Point 04

USI = [si_us0, cp_us0;
       si_us1, cp_us1;
       si_us2, cp_us2;
       si_us3, cp_us3;
       si_us4, cp_us4];

cp_us5 = ppval(pchip(USI(:,1),USI(:,2)),si_us5);                           % Upper Surface TE Cp Point 05 

si_ls0 = 0;                                                                % Lower Surface Arc Point 00
si_ls1 = si_ls0+sqrt((xi(4)-xi(5))^2+(zi(4)-zi(5))^2);                     % Lower Surface Arc Point 01
si_ls2 = si_ls1+sqrt((xi(3)-xi(4))^2+(zi(3)-zi(4))^2);                     % Lower Surface Arc Point 02
si_ls3 = si_ls2+sqrt((xi(2)-xi(3))^2+(zi(2)-zi(3))^2);                     % Lower Surface Arc Point 03
si_ls4 = si_ls3+sqrt((xi(1)-xi(2))^2+(zi(1)-zi(2))^2);                     % Lower Surface Arc Point 04
si_ls5 = si_ls4+sqrt((xj(1)-xi(1))^2+(zj(1)-zi(1))^2);                     % Lower Surface Arc Point 05

cp_ls0 = cpi(5);                                                           % Lower Surface Cp Point 00
cp_ls1 = cpi(4);                                                           % Lower Surface Cp Point 01
cp_ls2 = cpi(3);                                                           % Lower Surface Cp Point 02
cp_ls3 = cpi(2);                                                           % Lower Surface Cp Point 03
cp_ls4 = cpi(1);                                                           % Lower Surface Cp Point 04             

LSI = [si_ls0, cp_ls0;
       si_ls1, cp_ls1;
       si_ls2, cp_ls2;
       si_ls3, cp_ls3;
       si_ls4, cp_ls4];
   
cp_ls5 = ppval(pchip(LSI(:,1),LSI(:,2)),si_ls5);                           % Lower Surface TE Cp Point 05

arc_int_shift = sqrt((xj(1)-xi(1))^2+(zj(1)-zi(1))^2);                     % Initial Arc Shift 
arc_end_shift = sqrt((xj(end)-xi(end))^2+(zj(end)-zi(end))^2);             % Initial Arc Shift

% Note: As the si vector starts at 0 and sj begins a distance prior to this 
% all components of si need to be shifted to accomodate for the new arc 
% length incurred at the beginning of the arc. 'Arc_end_shift' is simply 
% added onto the final value of the vector so that it extends to the true 
% TE of the airfoil   

AF_XU = [0,                                        cp_ls5;
         si+arc_int_shift,                         cpi;
         si(end)+arc_int_shift+arc_end_shift,      cp_us5];

cpj = ppval(pchip(AF_XU(:,1),AF_XU(:,2)),sj);                              % Extrapolate the Cp for all panel end points

% Note: These additional lines of code force the coefficient of pressure to
% reach a value of 1 when the located at the stagnation point for 
% incompressible inviscid flow. This is discussed/proved in Moran P194. 

sj_pre = find(le(sj,SP(3)));
sj_pos = find(gt(sj,SP(3)));

sjsp = [sj(sj_pre); SP(3); sj(sj_pos)];
cjsp = [cpj(sj_pre); SP(5); cpj(sj_pos)];

cpj = ppval(pchip(sjsp,cjsp),sj);
end