function [outj] = fun_extrap(xi,zi,si,ini,xj,zj,sj)

% Note: The input value currently cannot be modelled against the arc length 
% because it is currently referenced against the xi components, which is the
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

in_us0 = ini(end-4);                                                       % Upper Surface Velocity Point 00                      
in_us1 = ini(end-3);                                                       % Upper Surface Velocity Point 01
in_us2 = ini(end-2);                                                       % Upper Surface Velocity Point 02
in_us3 = ini(end-1);                                                       % Upper Surface Velocity Point 03
in_us4 = ini(end);                                                         % Upper Surface Velocity Point 04

USI = [si_us0, in_us0;
       si_us1, in_us1;
       si_us2, in_us2;
       si_us3, in_us3;
       si_us4, in_us4];

in_us5 = ppval(pchip(USI(:,1),USI(:,2)),si_us5);                          % Upper Surface TE Velocity Point 05 

si_ls0 = 0;                                                                % Lower Surface Arc Point 00
si_ls1 = si_ls0+sqrt((xi(4)-xi(5))^2+(zi(4)-zi(5))^2);                     % Lower Surface Arc Point 01
si_ls2 = si_ls1+sqrt((xi(3)-xi(4))^2+(zi(3)-zi(4))^2);                     % Lower Surface Arc Point 02
si_ls3 = si_ls2+sqrt((xi(2)-xi(3))^2+(zi(2)-zi(3))^2);                     % Lower Surface Arc Point 03
si_ls4 = si_ls3+sqrt((xi(1)-xi(2))^2+(zi(1)-zi(2))^2);                     % Lower Surface Arc Point 04
si_ls5 = si_ls4+sqrt((xj(1)-xi(1))^2+(zj(1)-zi(1))^2);                     % Lower Surface Arc Point 05

in_ls0 = ini(5);                                                           % Lower Surface Velocity Point 00
in_ls1 = ini(4);                                                           % Lower Surface Velocity Point 01
in_ls2 = ini(3);                                                           % Lower Surface Velocity Point 02
in_ls3 = ini(2);                                                           % Lower Surface Velocity Point 03
in_ls4 = ini(1);                                                           % Lower Surface Velocity Point 04             

LSI = [si_ls0, in_ls0;
       si_ls1, in_ls1;
       si_ls2, in_ls2;
       si_ls3, in_ls3;
       si_ls4, in_ls4];
   
in_ls5 = ppval(pchip(LSI(:,1),LSI(:,2)),si_ls5);                       % Lower Surface TE Velocity Point 05

arc_int_shift = sqrt((xj(1)-xi(1))^2+(zj(1)-zi(1))^2);                     % Initial Arc Shift 
arc_end_shift = sqrt((xj(end)-xi(end))^2+(zj(end)-zi(end))^2);             % Initial Arc Shift

% Note: As the si vector starts at 0 and sj begins a distance prior to this 
% all components of si need to be shifted to accomodate for the new arc 
% length incurred at the beginning of the arc. 'Arc_end_shift' is simply 
% added onto the final value of the vector so that it extends to the true 
% TE of the airfoil   
% size(ini)
% size(si)
AF_XU = [0,                                        in_ls5;
         si+arc_int_shift,                         ini;
         si(end)+arc_int_shift+arc_end_shift,      in_us5];

outj = ppval(pchip(AF_XU(:,1),AF_XU(:,2)),sj);                              % Extrapolate the velocity for all panel end points

end
