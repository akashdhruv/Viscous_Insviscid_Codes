function [x0] = BL_initial_estimate(xi,zi,si,xj,zj,sj,Na,Nw,const)
  
loadfile = 'aero_variables.mat';
load(loadfile,'Uinf','Vinf')

si = si';
sj = sj';

% XXI = size(xi)
% Note: xi,zi,xj,zj are referenceing the coordinates on the airfoil only.

[uei,cpi,mui] = dirichlet_doublet_source(xi,zi,xj,zj,Vinf,Uinf,Na);

%         figure;
%         plot(xi,uei,'-k')
         
%         figure;
%         plot(xi,cpi,'-r')         

[uej] = fun_extrap(xi,zi,si,uei,xj,zj,sj);                                 % Velocity Profile Extrapolation
 
[aej] = vel_grad(sj,uej);

[SP] = stagnation_point(xj, zj, sj, uej);                                  % Stagnation Point Location

[cpj] = CP_extrap(xi, zi, si, cpi, xj, zj, sj, SP);                        % Coefficient of Pressure Extrapolation     
        
[dstarj, thetaj,~,~,~,~] = Coefficient_Of_Friction(const, xj, zj, sj, uej, aej, SP);

thetajw = (thetaj(1)+thetaj(end))*ones(Nw-1,1);                            % Initial Momentum Thickness

thetaj = [thetaj; thetajw];

[muj] = fun_extrap(xi,zi,si,mui,xj,zj,sj);                                 % Doublet Strength Across Airfoil Surface
        
 mujw = ((muj(1)+muj(end))/2)*ones(Nw-1,1);                                  % Doublet Strength Along Wake

 muj =  [muj; mujw];                                                       % Initial Doublet Strength
        
[mdefj] = fun_mdef(uej,dstarj,xj);                                         % Mass Defect Across Airfoil Surface
        
 mdefw= ((mdefj(1)+mdefj(end))/2)*ones(Nw+1,1);                              % Mass Defect Across Wake 
              
 mdef = [mdefj; mdefw];                                                    % Initial Mass Defect
 
 size(mdef)
 size(muj)
 size(thetaj)
 x0 = [-mdef, muj, thetaj];
 