function [cdHS] = calc_cdHS(Hi)

in = length(Hi);
cdHS = zeros(in,1);

for i = 1:in
    if lt(Hi(i),4)
        cdHS(i) = 0.207+0.00205*(4-Hi(i))^5.5;
    else
        cdHS(i) = 0.207-0.003*(Hi(i)-4)^2;
    end
end