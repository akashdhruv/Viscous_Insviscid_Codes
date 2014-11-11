function [avgTHT] = calc_avgTHT(theta)

in = length(theta);
avgTHT = zeros(in-1,1);

for i = 1:in-1
        avgTHT(i) = 0.5*(theta(i+1)+theta(i));
end
end
