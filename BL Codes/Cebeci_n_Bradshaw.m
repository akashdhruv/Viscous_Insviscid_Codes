function [S, H] = Cebeci_n_Bradshaw(lambda)

% Written Christopher Blower

% Defines the Cebeci and Bradshaw correlation formulas that derives
% Thwaites's graphs for Laminar Boundary Layer

lambda(2) = (lambda(1)+lambda(3))/2;
for i = 1:length(lambda)
if gt(lambda(i),0)                                                         % lt(lambda(i),0.1) & gt(lambda(i),0)
    S(i) = 0.22+(1.57*lambda(i))-(1.8*lambda(i)^2);
    H(i) = 2.61-(3.75*lambda(i))+(5.24*lambda(i)^2);
    
elseif lambda(i) == 0
    S(i) = 0.22;
    H(i) = 2.61;
    
elseif lt(lambda(i),0) & gt(lambda(i),-0.1)
    S(i) = 0.22+(1.402*lambda(i))+((0.018*lambda(i))/(lambda(i)+0.107));
    H(i) = 2.088+(0.0731/(lambda(i)+0.14));
else
    S(i) = 0;
    H(i) = 0;
end

S = S';
H = H';

end
end

    