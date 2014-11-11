function [cf] = laminar_cf(dstar, theta, ue, s)

% cf = laminar coefficient of friciton
% dstar = Displacement Thickness
% theta = Momentum Thickness
% ue = Tangential Velocity
% s = arc length from stagnation point

du_ds = zeros(1,length(s));
dt_ds = zeros(1,length(s));

for i = 1:length(s)
        if i == 1
            ds = s(2)-s(1);
            du = ue(2)-ue(1);
            dt = theta(2)-theta(1);
            du_ds(:,i) = du/ds;                                            % Change in Velocity wrt to Arc Length
            dt_ds(:,i) = dt/ds;                                            % Change in Momentum Thickness wrt to Arc Length
        elseif i == length(s)
            ds = s(end-1)-s(end);
            du = ue(end-1)-ue(end);
            dt = theta(end-1)-theta(end);
            du_ds(:,i) = du/ds;
            dt_ds(:,i) = dt/ds;
        else
            ds = s(i)-s(i-1);
            du = ue(i)-ue(i-1);
            dt = theta(i)-theta(i-1);
            du_ds(:,i) = du/ds;
            dt_ds(:,i) = dt/ds;
        end
end

in1 = (theta./ue)';
in2 = (2+(dstar./theta))';
in3 = (in1.*in2.*du_ds);

cf = 2*(dt_ds+in3);                                                        % von Karman momentum integral equation

end