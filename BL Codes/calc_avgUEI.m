function [avgUEI] = calc_avgUEI(Uei)

in = length(Uei);
avgUEI = zeros(in-1,1);

for i = 1:in-1
    avgUEI(i) = 0.5*(Uei(i+1)+Uei(i));
end
end