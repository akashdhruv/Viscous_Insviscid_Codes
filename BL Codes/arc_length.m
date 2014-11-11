function [s] = arc_length(x,z)

s(1) = 0;

for i = 2:length(x)
    s(i) = s(i-1)+(sqrt((x(i)-x(i-1))^2+(z(i)-z(i-1))^2));
end
