function [avgCF2] = calc_avgCF2(cf_2)

in = length(cf_2);
avgCF2 = zeros(in-1,1);

for i = in-1
    avgCF2(i) = 0.5*(cf_2(i+1)+cf_2(i));
end
end