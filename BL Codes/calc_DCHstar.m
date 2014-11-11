function [dcHS] = calc_DCHstar(Hi,Na)

in = length(Hi);
dcHS = zeros(in,1);

for i = 1:in
    if i >Na
        dcHS(i) = 1.52*((Hi(i)-1)^2)/(3+Hi(i)^3);
    else
        if lt(Hi(i),4)
            dcHS(i) = 0.207+0.00205*(4-Hi(i))^5.5;
        else
            dcHS(i) = 0.207-0.003*(Hi(i)-4)^2;
        end
    end
end
end