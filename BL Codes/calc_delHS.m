function [delHS] = calc_delHS(Hstar)

in = length(Hstar);
delHS = zeros(in-1,1);

for i = 1:in-1
    delHS(i) = Hstar(i+1)-Hstar(i); 
end
end