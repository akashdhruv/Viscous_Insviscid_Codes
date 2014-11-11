function [avgHS] = calc_avgHS(Hstar)

in = length(Hstar);
avgHS = zeros(in-1,1);

for i = 1:in-1
    avgHS(i) = 0.5*(Hstar(i+1)-Hstar(i));
end

