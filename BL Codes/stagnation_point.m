function [stag] = stagnation_point(xj,zj,sj,ue)
%% Velocity Profiling & Stagnation Point Identification
[ue_W ue_H] = size(ue);                                                    % Matrix Dimension Check - Confirm [001, 119] 
[xj_W xj_H] = size(xj);                                                    % Matrix Dimension Check - Confirm [001, 119]
[zj_W zj_H] = size(zj);                                                    % Matrix Dimension Check - Confirm [001, 119]

%         figure;
%         plot(x_DTA,ue_DTA,'-r')
%         title('Velocity Potential Distribution')
%         xlabel('Chordwise Location [%]')
%         ylabel('Velocity [m/s]')
%         axis square
%         axis([-0.1 1.1 -40 40])

%% Upper & Lower Surface Arc length

x_pp = pchip(sj,xj);
z_pp = pchip(sj,zj);
u_pp = pchip(sj,ue);

%         figure;
%         plot(sj,ppval(x_pp,sj),'-k')
%         axis equal
%         title('X-Coordinate wrt Arc Length') 
%         xlabel('Arc Length')
%         ylabel('X-Position')
%         
%         figure;
%         plot(sj,ppval(z_pp,sj),'-k')
%         axis equal
%         title('Z-Coordinate wrt Arc Length') 
%         xlabel('Arc Length')
%         ylabel('Z-Position')
%         
%         figure;
%         plot(sj,ppval(u_pp,sj),'-k')
%         title('Tangential Velocity Profile wrt Arc Length') 
%         xlabel('Arc Length')
%         ylabel('Velocity')

LL = floor(length(sj)*(1/4));                                              % By multiplying by 1/4 & 3/4 respectively this removes the
UL = ceil(length(sj)*(3/4));                                               % Trailing edge of the airfoil, as we want to assess the airfoil 
                                                                           % about the LE. If we have a stagnation point at the TE => Huge Error!

sj_sp = fzero(@fn1,[sj(LL),sj(UL)]);                                       % Arc Length Location of the Stagnation Point

    function u = fn1(sj)
        u = ppval(u_pp,sj);
    end

xj_sp = ppval(x_pp,sj_sp);
zj_sp = ppval(z_pp,sj_sp);
uj_sp = 0;
cj_sp = 1;

stag = [xj_sp, zj_sp, sj_sp, uj_sp, cj_sp];                                % Stagnation Point Location
end
