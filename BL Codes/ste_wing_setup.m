function [xi,zi,si,xj,zj,sj,n_com,t_com,Nw,Naw,del_BSr,xia,zia,sia,xja,zja,sja] = ste_wing_setup(NACA,c,Na,Nw)

%% Airfoil Coordinates
Ns = Na/2;                                                                 % Number of Panels per Surface

[xu,zu,xl,zl] = naca4m(NACA,c,Ns);                                         % NACA Airfoil Generator

% figure;
% plot(xu,zu,'-b',...
%      xl,zl,'-r')
% axis equal
% axis([0 1.4 -0.5 0.5])

xja = [xl; flipud(xu(1:end-1))];                                           % Body - X-coord. of End Points
zja = [zl; flipud(zu(1:end-1))];                                           % Body - Z-coord. of End Points
sja = arc_length(xja,zja);                                                 % Body - Arc Length

xjw = (linspace(1,11,Nw))';                                               % Wake - X-coord. of End Points (Horizontal Wake Employed)
zjw = zeros(Nw,1);                                                         % Wake - Z-coord. of End Points

xj = [xja; xjw(2:end)];
zj = [zja; zjw(2:end)];

sj = arc_length(xj,zj);

Nw = Nw-1;                      
Naw= Na+Nw;                                                                % Total Number of Panels for Airfoil and Wake

xia= (xja(2:end)+xja(1:end-1))/2;                                          % Body - X-coord. of Collocation Point
zia= (zja(2:end)+zja(1:end-1))/2;                                          % Body - Z-coord. of Collocation Point
sia= arc_length(xia,zia);

xi = (xj(2:end)+xj(1:end-1))/2;                                            % Body & Wake - X-coord. of Collocation Point
zi = (zj(2:end)+zj(1:end-1))/2;                                            % Body & Wake - Z-coord. of Collocation Point 

si = arc_length(xi,zi);

for n = 1:Naw
	[n_com(:,n),t_com(:,n)] = n_t_vectors(xj(n),zj(n),xj(n+1),zj(n+1));
end

% figure;
% plot(xj,zj,'b*-',xi,zi,'rx');
% axis equal;

del_BSr = 0;

end