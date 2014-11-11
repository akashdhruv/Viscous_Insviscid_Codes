function [uei,cpi,mui] = dirichlet_doublet_source(xi,zi,xj,zj,Vinf,Uinf,N)

n_vec = zeros(2,N);
t_vec = zeros(2,N);

for n = 1:N
	[n_vec(:,n),t_vec(:,n)] = n_t_vectors(xj(n),zj(n),xj(n+1),zj(n+1));
end

%% Calculate Influence Matrices
n = zeros(2,N);
A = zeros(N,N);
B = zeros(N,N); 
C = zeros(N+1,N+1);    

wake_x = 10^7;                                                             % X Length of Wake Geometry
wake_z = 0;                                              % Z Length of Wake Geometry

for i = 1:N
   for j = 1:N
        if i == j 
            B(i,j) = PHICS(1,xi(i),zi(i),xj(j),zj(j),xj(j+1),zj(j+1));
            C(i,j) = 0.5; 
        else
            B(i,j) = PHICS(1,xi(i),zi(i),xj(j),zj(j),xj(j+1),zj(j+1));
            C(i,j) = PHICD(1,xi(i),zi(i),xj(j),zj(j),xj(j+1),zj(j+1));
        end
   end
   C(i,N+1) = PHICD(1,xi(i),zi(i),xj(1),zj(1),wake_x,wake_z+zj(1));
end

%% Trailing Edge Closure Conditions
C(N+1,1)  = 1;                                                             
C(N+1,N)  = -1;
C(N+1,N+1)= 1;                                                              

%% Convert to the A Matrix 
for i = 1:N
    for j = 1:N
        if eq(j,1)
            A(i,1) = C(i,1)-C(i,N+1);
        elseif eq(j,N)
            A(i,N) = C(i,N)+C(i,N+1);
        else
            A(i,j) = C(i,j);
        end
    end
end

RHS_DO = zeros(N+1,1);
SIG_SD = zeros(N,1);

for i = 1:N
    SIG_SD(i,1) = [Uinf Vinf]*n_vec(:,i);
    RHS_DO(i,1)   = -[xi(i) zi(i)]*[Uinf;Vinf];                            % Right Hand Side Vector - Doublet Only
end

RHS_SD =(-B)*SIG_SD;                                                       % Right Hand Side Vector - Source Doublet

MU_DO = C\RHS_DO;                                                          % Doublet Strength - Doublet Method Only(mu)
MU_SD = A\RHS_SD;                                                          % Doublet Strength - Source-Doublet Method Only(mu)

for i = 1:N
    v_inf_t = [Uinf Vinf]*t_vec(:,i);
        if i == 1
            RR = sqrt((xi(2)-xi(1))^2+(zi(2)-zi(1))^2);
            v_loc_t_DO = (MU_DO(2)-MU_DO(1))/RR;                           % Local Tangential Velocity - Doublet Method Only
            v_loc_t_SD = (MU_SD(2)-MU_SD(1))/RR;                           % Local Tangential Velocity - Source-Doublet Method Only
        elseif i == N
            RR = sqrt((xi(N)-xi(N-1))^2+(zi(N)-zi(N-1))^2);
            v_loc_t_DO = (MU_DO(N)-MU_DO(N-1))/RR;                         % Local Tangential Velocity - Doublet Method Only
            v_loc_t_SD = (MU_SD(N)-MU_SD(N-1))/RR;                         % Local Tangential Velocity - Source-Doublet Method Only
        else
            RR = sqrt((xi(i+1)-xi(i-1))^2+(zi(i+1)-zi(i-1))^2);
            v_loc_t_DO = (MU_DO(i+1)-MU_DO(i-1))/RR;                       % Local Tangential Velocity - Doublet Method Only
            v_loc_t_SD = ((MU_SD(i+1)-MU_SD(i-1))/RR);                     % Local Tangential Velocity - Source-Doublet Method Only
        end
    
    Ue_SDS(i) = v_loc_t_SD + v_inf_t;
    Cp_CSD(i) = 1 - ((v_loc_t_DO)^2/Uinf^2);                               % Coefficient of Pressure Doublet Only
    Cp_SDS(i) = 1 - ((v_loc_t_SD + v_inf_t)^2/Uinf^2);                     % Coefficient of Pressure Source-Doublet
end
mui = MU_SD;
uei = Ue_SDS';
cpi = Cp_SDS';

% UUU = size(uei)
% CCC = size(cpi)
% MMM = size(mui)
ind = length(xi)/2;

%% Plot Pressure Profiles

% figure;
% plot(xi,Cp_CSD,'-b',xi,Cp_SDS,'-g');
% set(gca,'YDir','reverse');
% title('C_p for Doublet and Source-Doublet Panel Method Code across a NACA4412')
% ylabel('Coefficient of Pressure, C_p [-]');
% xlabel('Normalized Chord Length, c [-]');
% legend('Doublet Only','Source-Doublet');
% axis([0 1 -2 2])
end
