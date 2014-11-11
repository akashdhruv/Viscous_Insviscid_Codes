function [cf_2] = calc_cf(Hi,Na)

in = length(Hi);
cf_2 = zeros(in,1);

for i = 1:in
    if i >Na
        cf_2(i) = 0;
    else
        if lt(Hi,7.4)
            cf_2(i) = -0.067+0.01977*(((7.4-Hi(i))^2)/(Hi(i)-1));
        else 
            cf_2(i) = -0.067+0.022*(1-(1.4/(Hi(i)-6)))^2;
        end
    end
end
end