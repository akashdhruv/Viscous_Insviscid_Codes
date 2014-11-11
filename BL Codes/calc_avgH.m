function [H] = calc_avgH(Hi)

in = length(Hi);
H = zeros(in-1,1);

for i = 1:in-1;
    H(i) = 0.5*(Hi(i+1)+Hi(i));
end
end