function [delX] = calc_delX(xi)

in = length(xi);
delX = zeros(in-1,1);

for i = 1:in-1
    delX(i) = xi(i+1)-xi(i);
end
end
 