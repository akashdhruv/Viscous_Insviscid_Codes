function [delUEI] = calc_delUEI(Uei)

in = length(Uei);
delUEI = zeros(in-1,1);

for i = 1:in-1
    delUEI(i) = Uei(i+1)-Uei(i);
end
end