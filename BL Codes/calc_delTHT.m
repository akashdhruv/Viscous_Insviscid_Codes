function [delTHT] = calc_delTHT(theta)

in = length(theta);
delTHT = zeros(in-1,1);

for i = 1:in-1
    delTHT(i) = theta(i+1)-theta(i);
end
end