function [avgDCH] = calc_avgDCH(DC_Hstar)

in = length(DC_Hstar);
avgDCH = zeros(in-1,1);

for i = 1:in-1
    avgDCH(i) = 0.5*(DC_Hstar(i+1)+DC_Hstar(i));
end
end